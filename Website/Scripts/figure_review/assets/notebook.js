import { renderMarkdown } from './markdown.js';

const $ = (id) => document.getElementById(id);
const pageLinks = [...document.querySelectorAll('[data-page]')];
const paper = $('paper');
const figureUrl = (path) => `/figure?${new URLSearchParams({ figure: path })}`;
let current = null;
let figures = [];
let libraryMessage = 'Loading figures…';
let loading = false;
let preview = false;
let insertionIndex = null;
let drag = null;
let timer;

function element(tag, className, label) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (label) node.textContent = label;
  return node;
}

function control(label, action, className = 'button') {
  const button = element('button', className, label);
  button.type = 'button';
  button.addEventListener('click', action);
  return button;
}

async function request(path, body) {
  const response = await fetch(path, body && {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body),
  });
  const result = await response.json();
  if (!response.ok) throw Object.assign(new Error(result.error), { status: response.status });
  return result;
}

function status(message, isError = false) {
  $('save-status').textContent = message;
  $('save-status').classList.toggle('error', isError);
}

function showAccess() {
  $('access').hidden = false;
  $('workspace').hidden = true;
  $('code').focus();
}

function updateControls() {
  $('save').disabled = !current?.dirty || Boolean(current?.saving) || Boolean(current?.conflict) || loading;
  $('preview').disabled = !current || loading;
  $('reload').hidden = !current?.conflict && (Boolean(current) || loading);
  for (const button of pageLinks) button.disabled = loading;
  for (const button of $('figures').querySelectorAll('button')) button.disabled = !current || loading || preview;
}

function changed() {
  current.dirty = true;
  current.change += 1;
  if (!current.conflict) {
    status('Unsaved changes');
    clearTimeout(timer);
    timer = setTimeout(() => void save(current), 700);
  }
  updateControls();
}

async function save(page) {
  clearTimeout(timer);
  if (!page || page.conflict) return false;
  if (page.saving) {
    if (!await page.saving) return false;
    return save(page);
  }
  if (!page.dirty) return true;
  const version = page.change;
  status('Saving…');
  page.saving = (async () => {
    try {
      const { revision } = await request('/api/notebook', {
        page: page.id, revision: page.revision, blocks: page.blocks,
      });
      page.revision = revision;
      page.dirty = version !== page.change;
      status(page.dirty ? 'Unsaved changes' : 'Saved');
      return true;
    } catch (error) {
      page.conflict = error.status === 409;
      status(error.message, true);
      if (error.status === 401) showAccess();
      return false;
    } finally {
      page.saving = null;
      updateControls();
    }
  })();
  updateControls();
  return page.saving;
}

async function openPage(id, discard = false) {
  if (loading || (!discard && current?.id === id)) return;
  loading = true;
  updateControls();
  try {
    while (!discard && current?.dirty) {
      if (!await save(current)) return;
    }
    const { revision, blocks } = await request(`/api/notebook?${new URLSearchParams({ page: id })}`);
    current = { id, revision, blocks, change: 0, dirty: false, saving: null, conflict: false };
    insertionIndex = null;
    const selected = pageLinks.find((button) => button.dataset.page === id);
    selected.closest('details').open = true;
    for (const button of pageLinks) button.removeAttribute('aria-current');
    selected.setAttribute('aria-current', 'page');
    $('page-title').textContent = selected.textContent;
    $('page-group').textContent = selected.closest('details').querySelector('summary').textContent;
    history.replaceState(null, '', `#${id}`);
    $('access').hidden = true;
    $('workspace').hidden = false;
    status('');
    renderPage();
  } catch (error) {
    if (error.status === 401) showAccess();
    else {
      $('workspace').hidden = false;
      $('access').hidden = true;
      status(error.message, true);
    }
  } finally {
    loading = false;
    updateControls();
  }
}

function finishDrag() {
  drag = null;
  paper.classList.remove('dragging');
  paper.querySelector('.drop-active')?.classList.remove('drop-active');
}

function draggable(node, data) {
  const targetAt = (event) => document.elementFromPoint(event.clientX, event.clientY)?.closest('.drop-line');
  node.addEventListener('pointerdown', (event) => {
    if (!current || loading || preview || drag || event.button !== 0 || event.target.closest('button')) return;
    event.preventDefault();
    node.setPointerCapture(event.pointerId);
    drag = { ...data, pointer: event.pointerId };
  });
  node.addEventListener('pointermove', (event) => {
    if (drag?.pointer !== event.pointerId) return;
    paper.classList.add('dragging');
    paper.querySelector('.drop-active')?.classList.remove('drop-active');
    targetAt(event)?.classList.add('drop-active');
  });
  node.addEventListener('pointerup', (event) => {
    if (drag?.pointer !== event.pointerId) return;
    const item = drag;
    const index = targetAt(event)?.dataset.index;
    finishDrag();
    if (index === undefined || loading) return;
    if (item.figure) insertFigure(item.figure, Number(index));
    else moveBlock(item.index, Number(index));
  });
  node.addEventListener('lostpointercapture', (event) => {
    if (drag?.pointer === event.pointerId) finishDrag();
  });
}

function renderLibrary() {
  const filter = $('search').value.trim().toLowerCase();
  const matches = figures.filter((path) => path.toLowerCase().includes(filter));
  $('figures').replaceChildren();
  if (!matches.length) {
    $('figures').append(element('p', 'muted', libraryMessage || 'No matching figures.'));
    return;
  }
  for (const figure of matches) {
    const row = element('div', 'drive-row');
    const path = element('span', 'drive-path', figure);
    const insert = control('Insert', () => insertFigure(figure, insertionIndex ?? current.blocks.length));
    insert.setAttribute('aria-label', `Insert ${figure}`);
    draggable(row, { figure });
    row.append(path, insert);
    $('figures').append(row);
  }
  updateControls();
}

async function refreshLibrary() {
  $('refresh').disabled = true;
  try {
    ({ figures } = await request('/api/figures'));
    libraryMessage = figures.length ? '' : 'No figures yet. Add files in Drive and refresh.';
  } catch (error) {
    if (error.status === 401) showAccess();
    figures = [];
    libraryMessage = error.message;
  } finally {
    $('refresh').disabled = false;
    renderLibrary();
  }
}

function insertBlock(block, index) {
  if (!current || loading || preview) return;
  if (current.blocks.length >= 200) {
    status('This page can hold up to 200 blocks.', true);
    return;
  }
  current.blocks.splice(index, 0, block);
  insertionIndex = index + 1;
  renderPage();
  changed();
  paper.querySelectorAll('.document-block')[index].querySelector('textarea')?.focus();
}

function insertFigure(figure, index) {
  insertBlock({ kind: 'figure', figure, layout: 'full', markdown: '' }, index);
}

function moveBlock(from, to) {
  if (loading || from === to || from + 1 === to) return;
  const [block] = current.blocks.splice(from, 1);
  const target = from < to ? to - 1 : to;
  current.blocks.splice(target, 0, block);
  insertionIndex = target + 1;
  renderPage();
  changed();
}

function dropLine(index) {
  const line = element('div', 'drop-line');
  line.dataset.index = index;
  const add = control('Add text', () => insertBlock({ kind: 'text', markdown: '' }, index));
  add.setAttribute('aria-label', `Add text at position ${index + 1}`);
  line.append(add);
  return line;
}

function textContent(block, index) {
  if (preview) return renderMarkdown(block.markdown);
  const field = element('textarea', 'block-text');
  field.rows = block.kind === 'figure' ? 3 : 4;
  field.maxLength = 20_000;
  field.placeholder = block.kind === 'figure' ? 'Write about this figure in Markdown…' : 'Write in Markdown…';
  field.setAttribute('aria-label', block.kind === 'figure' ? 'Figure text' : 'Markdown text');
  field.value = block.markdown;
  field.addEventListener('focus', () => { insertionIndex = index + 1; });
  field.addEventListener('input', () => {
    block.markdown = field.value;
    changed();
  });
  return field;
}

function figureContent(block, text) {
  const wrapper = element('div', `figure-content ${block.layout}`);
  const media = element('figure', 'media-box');
  const name = block.figure.split('/').at(-1);
  const url = figureUrl(block.figure);
  const pdf = name.toLowerCase().endsWith('.pdf');
  const image = element(pdf ? 'iframe' : 'img');
  if (pdf) image.title = name;
  else {
    image.alt = name;
    image.loading = 'lazy';
    image.addEventListener('error', () => image.replaceWith(element('p', 'error', 'Figure unavailable in Results_Main.')));
  }
  image.src = url;
  const caption = element('figcaption', 'muted');
  const source = element('a', '', 'Open original');
  source.href = url;
  source.target = '_blank';
  source.rel = 'noopener';
  caption.append(element('span', '', name), source);
  media.append(image, caption);
  wrapper.append(media);
  if (!preview || block.markdown.trim()) wrapper.append(text);
  return wrapper;
}

function blockTools(block, index, content) {
  const tools = element('div', 'block-tools');
  if (block.kind === 'figure') {
    const handle = element('span', 'drag-handle', 'Drag to move');
    draggable(handle, { index });
    const layout = element('select');
    layout.setAttribute('aria-label', 'Figure layout');
    for (const [value, label] of [['full', 'Text below'], ['left', 'Figure left'], ['right', 'Figure right']]) {
      const option = element('option', '', label);
      option.value = value;
      layout.append(option);
    }
    layout.value = block.layout;
    layout.addEventListener('change', () => {
      if (loading) return;
      content.classList.replace(block.layout, layout.value);
      block.layout = layout.value;
      changed();
    });
    tools.append(handle, layout);
  }
  const up = control('Move up', () => moveBlock(index, index - 1));
  const down = control('Move down', () => moveBlock(index, index + 2));
  up.disabled = index === 0;
  down.disabled = index === current.blocks.length - 1;
  tools.append(up, down, control('Remove', () => {
    if (loading || ((block.markdown.trim() || block.kind === 'figure') && !window.confirm('Remove this block?'))) return;
    current.blocks.splice(index, 1);
    insertionIndex = index;
    renderPage();
    changed();
  }));
  return tools;
}

function renderPage() {
  paper.replaceChildren();
  if (!current.blocks.length) {
    paper.append(element('p', 'empty-hint', preview ? 'This page is empty.' : 'Start with text or drag in a figure from Drive.'));
  }
  current.blocks.forEach((block, index) => {
    if (!preview) paper.append(dropLine(index));
    const article = element('section', preview ? 'preview-block' : 'document-block');
    const text = textContent(block, index);
    const content = block.kind === 'figure' ? figureContent(block, text) : text;
    if (!preview) article.append(blockTools(block, index, content));
    article.append(content);
    paper.append(article);
  });
  if (!preview) paper.append(dropLine(current.blocks.length));
}

async function start() {
  let hash = location.hash.slice(1);
  if (hash.startsWith('invite=')) {
    const code = hash.slice('invite='.length);
    hash = '';
    history.replaceState(null, '', location.pathname + location.search);
    try {
      await request('/api/login', { code });
    } catch (error) {
      showAccess();
      $('login-error').textContent = error.status === 401 ? 'This invite is invalid or has expired.' : error.message;
      return;
    }
  }
  const id = pageLinks.some((button) => button.dataset.page === hash) ? hash : 'results-1';
  await openPage(id);
  if (current) await refreshLibrary();
}

for (const button of pageLinks) button.addEventListener('click', () => void openPage(button.dataset.page));
$('preview').addEventListener('click', () => {
  preview = !preview;
  $('preview').setAttribute('aria-pressed', String(preview));
  $('preview').textContent = preview ? 'Edit Markdown' : 'Preview Markdown';
  paper.classList.toggle('preview', preview);
  renderPage();
  updateControls();
});
$('save').addEventListener('click', () => void save(current));
$('reload').addEventListener('click', () => {
  if (!current) void start();
  else if (window.confirm('Discard local changes and reload this page?')) void openPage(current.id, true);
});
$('refresh').addEventListener('click', refreshLibrary);
$('search').addEventListener('input', renderLibrary);
$('login').addEventListener('submit', async (event) => {
  event.preventDefault();
  $('login-error').textContent = '';
  try {
    await request('/api/login', { code: $('code').value });
    $('login').reset();
    if (current) {
      $('access').hidden = true;
      $('workspace').hidden = false;
      await save(current);
      await refreshLibrary();
    } else await start();
  } catch (error) {
    $('login-error').textContent = error.message;
  }
});
window.addEventListener('beforeunload', (event) => {
  if (!current?.dirty) return;
  event.preventDefault();
  event.returnValue = '';
});
document.addEventListener('keydown', (event) => {
  if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === 's' && current && !$('workspace').hidden) {
    event.preventDefault();
    void save(current);
  }
});

void start();
