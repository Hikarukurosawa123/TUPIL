import { renderMarkdown } from './markdown.js';
import { blockId, sourceAnchor, previewAnchor, markComments, revealAnchor } from './anchors.js';

const $ = (id) => document.getElementById(id);
const sections = [...document.querySelectorAll('[data-section]')];
const paper = $('paper');
const figureUrl = (path) => `/figure?${new URLSearchParams({ figure: path })}`;
const commentDrafts = new Map();
const commentTime = new Intl.DateTimeFormat(undefined, { dateStyle: 'medium', timeStyle: 'short' });
let current = null;
let pages = [];
let trash = [];
let selectedAnchor = null;
let createSection = null;
let organizing = false;
let pageLoadFailed = false;
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
  $('reload').hidden = !current?.conflict && (Boolean(current) || loading || !pageLoadFailed);
  $('rename').hidden = !current;
  $('rename').disabled = loading || organizing || Boolean(current?.renaming);
  $('move-page').hidden = !current;
  $('move-page').disabled = loading || organizing || Boolean(current?.renaming);
  $('delete-page').hidden = !current;
  $('delete-page').disabled = loading || organizing || Boolean(current?.renaming || current?.postingComment);
  $('comments-jump').disabled = !current || loading;
  for (const button of document.querySelectorAll('.page-link, .new-page, .restore-page')) {
    button.disabled = loading || organizing || Boolean(current?.renaming || current?.postingComment);
  }
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
  if (loading || organizing || current?.renaming || current?.postingComment || (!discard && current?.id === id)) return;
  loading = true;
  pageLoadFailed = false;
  updateControls();
  try {
    while (!discard && current?.dirty) {
      if (!await save(current)) return;
    }
    const { revision, blocks, pages: nextPages, trash: nextTrash } = await request(`/api/notebook?${new URLSearchParams({ page: id })}`);
    current = { id, revision, blocks, change: 0, dirty: false, saving: null, conflict: false, comments: null, commentRequest: 0 };
    insertionIndex = null;
    renderNavigation(nextPages, nextTrash);
    history.replaceState(null, '', `#${id}`);
    $('access').hidden = true;
    $('workspace').hidden = false;
    closeRename();
    closePageForm('create');
    closePageForm('move');
    status('');
    renderPage();
    $('page-comments').hidden = false;
    showCommentDraft();
    commentStatus('Loading comments…');
    renderComments(current);
    void loadComments(current);
  } catch (error) {
    if (error.status === 401) showAccess();
    else {
      $('workspace').hidden = false;
      $('access').hidden = true;
      pageLoadFailed = true;
      status(error.message, true);
    }
  } finally {
    loading = false;
    updateControls();
  }
}

function sectionLabel(section) {
  return sections.find((item) => item.dataset.section === section).querySelector('summary').textContent;
}

function renderNavigation(nextPages, nextTrash = trash) {
  pages = nextPages;
  trash = nextTrash;
  for (const section of sections) {
    const list = section.querySelector('.page-links');
    list.replaceChildren();
    for (const page of pages.filter((entry) => entry.section === section.dataset.section)) {
      const button = control(page.title, () => void openPage(page.id), 'page-link');
      if (current?.id === page.id) {
        button.setAttribute('aria-current', 'page');
        section.open = true;
        $('page-title').textContent = page.title;
        $('page-group').textContent = sectionLabel(page.section);
      }
      list.append(button);
    }
  }
  $('trash').hidden = !trash.length;
  $('trash-pages').replaceChildren();
  for (const page of trash) {
    const row = element('li');
    const restore = control('Restore', () => void changeTrash(page.id, 'restore'), 'button restore-page');
    restore.setAttribute('aria-label', `Restore ${page.title}`);
    row.append(element('span', '', page.title), restore);
    $('trash-pages').append(row);
  }
  updateControls();
}

function showEmptyWorkspace() {
  current = null;
  pageLoadFailed = false;
  clearSelection();
  closeRename();
  closePageForm('move');
  $('access').hidden = true;
  $('workspace').hidden = false;
  $('page-title').textContent = 'No pages';
  $('page-group').textContent = '';
  $('page-comments').hidden = true;
  paper.replaceChildren(element('p', 'empty-hint', 'Create a page from the sidebar or restore one from Trash.'));
  markComments(paper, []);
  history.replaceState(null, '', location.pathname + location.search);
  updateControls();
  $('reload').hidden = true;
}

async function changeTrash(id, action) {
  if (loading || organizing || current?.renaming || current?.postingComment) return;
  const deleting = action === 'delete';
  if (deleting && !window.confirm('Move this page to Trash? You can restore its content and comments.')) return;
  const section = pages.find((page) => page.id === id)?.section;
  organizing = true;
  updateControls();
  try {
    while (deleting && current?.dirty) {
      if (!await save(current)) return;
    }
    const { pages: remaining, trash: deleted } = await request('/api/notebook/trash', { page: id, action });
    if (deleting) {
      commentDrafts.delete(id);
      showEmptyWorkspace();
    }
    renderNavigation(remaining, deleted);
    organizing = false;
    if (deleting) {
      const next = remaining.find((page) => page.section === section) ?? remaining[0];
      if (next) await openPage(next.id);
      $('trash').open = true;
      status('Page moved to Trash.');
    } else {
      await openPage(id);
      if (current?.id === id) status('Page restored.');
    }
  } catch (error) {
    if (error.status === 401) showAccess();
    status(error.message, true);
  } finally {
    organizing = false;
    updateControls();
  }
}

function closePageForm(action) {
  $(`${action}-page-form`).hidden = true;
  $(`${action}-page-status`).textContent = '';
}

function showCreate(section) {
  closeRename();
  closePageForm('move');
  closePageForm('create');
  createSection = section.dataset.section;
  const form = $('create-page-form');
  section.after(form);
  $('create-page-label').textContent = `New page title in ${sectionLabel(createSection)}`;
  $('new-page-title').value = '';
  form.hidden = false;
  $('new-page-title').focus();
}

function showMove() {
  closeRename();
  closePageForm('create');
  closePageForm('move');
  const page = pages.find((item) => item.id === current.id);
  $('move-section').value = page.section;
  const inSection = pages.filter((item) => item.section === page.section);
  const previous = inSection[inSection.findIndex((item) => item.id === page.id) - 1]?.id ?? '';
  movePositions();
  $('move-position').value = previous;
  $('move-page-form').hidden = false;
  $('move-section').focus();
}

function movePositions() {
  const select = $('move-position');
  select.replaceChildren(new Option('First in section', ''));
  for (const page of pages.filter((item) => item.section === $('move-section').value && item.id !== current.id)) {
    select.add(new Option(`After ${page.title}`, page.id));
  }
  select.selectedIndex = select.options.length - 1;
}

async function organizePages(action, endpoint, body) {
  const form = $(`${action}-page-form`);
  organizing = true;
  updateControls();
  $(`${action}-page-status`).textContent = '';
  const inputs = form.querySelectorAll('input, select, button');
  for (const control of inputs) control.disabled = true;
  try {
    const { pages: updated, trash: deleted, id } = await request(endpoint, body);
    renderNavigation(updated, deleted);
    closePageForm(action);
    return id ?? current?.id;
  } catch (error) {
    if (error.status === 401) showAccess();
    $(`${action}-page-status`).textContent = error.message;
    return null;
  } finally {
    organizing = false;
    for (const control of inputs) control.disabled = false;
    updateControls();
  }
}

function closeRename() {
  $('rename-form').hidden = true;
  $('rename-status').textContent = '';
}

async function renamePage(event) {
  event.preventDefault();
  if (!current || loading || current.renaming) return;
  const page = current;
  page.renaming = true;
  $('rename-status').textContent = '';
  for (const id of ['title-input', 'rename-cancel', 'rename-save']) $(id).disabled = true;
  updateControls();
  try {
    const { title } = await request('/api/notebook/title', { page: page.id, title: $('title-input').value });
    renderNavigation(pages.map((entry) => entry.id === page.id ? { ...entry, title } : entry));
    closeRename();
  } catch (error) {
    if (error.status === 401) showAccess();
    $('rename-status').textContent = error.message;
  } finally {
    page.renaming = false;
    for (const id of ['title-input', 'rename-cancel', 'rename-save']) $(id).disabled = false;
    updateControls();
    if ($('rename-form').hidden) $('rename').focus();
  }
}

function commentStatus(message = '', error = false) {
  $('comment-status').textContent = message;
  $('comment-status').classList.toggle('error', error);
}

function offerSelection(anchor) {
  selectedAnchor = anchor && anchor.quote.length <= 2_000 ? anchor : null;
  $('selection-tools').hidden = !selectedAnchor;
}

function clearSelection() {
  offerSelection(null);
}

function showCommentDraft() {
  const draft = current && commentDrafts.get(current.id);
  $('comment-form').hidden = !draft;
  $('comment-help').hidden = Boolean(draft);
  $('comment-quote').textContent = draft?.anchor.quote ?? '';
  $('comment-text').value = draft?.text ?? '';
}

function commentOnSelection() {
  if (!selectedAnchor || !current || current.postingComment) return;
  const previous = commentDrafts.get(current.id);
  commentDrafts.set(current.id, { anchor: selectedAnchor, text: previous?.text ?? '' });
  clearSelection();
  showCommentDraft();
  $('comment-text').focus();
}

function showCommentSelection(anchor) {
  setPreview(anchor.mode === 'preview');
  if (!revealAnchor(paper, anchor)) commentStatus('The original selection has changed or been removed. Its quote is still shown below.', true);
}

function renderComments(page) {
  $('comments-list').replaceChildren();
  $('comments-empty').hidden = page.comments?.length !== 0;
  for (const comment of page.comments ?? []) {
    const item = element('li', 'comment-item');
    const meta = element('div', 'comment-meta');
    const time = element('time', 'muted', commentTime.format(new Date(comment.created)));
    time.dateTime = comment.created;
    meta.append(element('strong', '', comment.author), time);
    item.append(meta);
    if (comment.anchor) item.append(element('blockquote', 'comment-quote', comment.anchor.quote));
    item.append(element('p', '', comment.text));
    if (comment.anchor) item.append(control('Show selection', () => showCommentSelection(comment.anchor)));
    $('comments-list').append(item);
  }
  markComments(paper, page.comments ?? []);
}

async function loadComments(page, successMessage = '') {
  const requestId = ++page.commentRequest;
  $('comments-refresh').disabled = true;
  try {
    const { comments } = await request(`/api/notebook/comments?${new URLSearchParams({ page: page.id })}`);
    if (current !== page || requestId !== page.commentRequest) return;
    page.comments = comments;
    renderComments(page);
    commentStatus(successMessage);
  } catch (error) {
    if (current !== page || requestId !== page.commentRequest) return;
    if (error.status === 401) showAccess();
    commentStatus(successMessage ? `${successMessage} Could not refresh the list.` : error.message, true);
  } finally {
    if (current === page && requestId === page.commentRequest) $('comments-refresh').disabled = Boolean(page.postingComment);
  }
}

async function postComment(event) {
  event.preventDefault();
  if (!current || loading || organizing || current.postingComment) return;
  const page = current;
  const draft = commentDrafts.get(page.id);
  if (!draft) return;
  page.postingComment = true;
  const controls = ['comment-text', 'comment-author', 'comment-post', 'comment-cancel', 'comments-refresh'];
  for (const id of controls) $(id).disabled = true;
  updateControls();
  commentStatus('Posting comment…');
  try {
    while (page.dirty) {
      if (!await save(page)) {
        commentStatus('Save the page before posting a comment on this selection.', true);
        return;
      }
    }
    const { comment } = await request('/api/notebook/comments', {
      page: page.id, author: $('comment-author').value, text: draft.text, anchor: draft.anchor,
    });
    commentDrafts.delete(page.id);
    showCommentDraft();
    page.comments = [comment, ...(page.comments ?? []).filter((item) => item.id !== comment.id)];
    renderComments(page);
    await loadComments(page, 'Comment posted.');
  } catch (error) {
    if (error.status === 401) showAccess();
    commentStatus(error.message, true);
  } finally {
    page.postingComment = false;
    for (const id of controls) $(id).disabled = false;
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
  current.blocks.splice(index, 0, { id: blockId(), ...block });
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
  if (preview) {
    const formatted = renderMarkdown(block.markdown);
    formatted.dataset.block = block.id;
    return formatted;
  }
  const field = element('textarea', 'block-text');
  field.rows = block.kind === 'figure' ? 3 : 4;
  field.maxLength = 20_000;
  field.placeholder = block.kind === 'figure' ? 'Write about this figure in Markdown…' : 'Write in Markdown…';
  field.setAttribute('aria-label', block.kind === 'figure' ? 'Figure text' : 'Markdown text');
  field.value = block.markdown;
  field.dataset.block = block.id;
  field.addEventListener('select', () => offerSelection(sourceAnchor(field)));
  field.addEventListener('pointerup', () => offerSelection(sourceAnchor(field)));
  field.addEventListener('focus', () => { insertionIndex = index + 1; });
  field.addEventListener('input', () => {
    clearSelection();
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
  clearSelection();
  paper.replaceChildren();
  if (!current.blocks.length) {
    paper.append(element('p', 'empty-hint', preview ? 'This page is empty.' : 'Start with text or drag in a figure from Drive.'));
  }
  current.blocks.forEach((block, index) => {
    if (preview && block.kind === 'text' && !block.markdown.trim()) return;
    if (!preview) paper.append(dropLine(index));
    const article = element('section', preview ? 'preview-block' : 'document-block');
    const text = textContent(block, index);
    const content = block.kind === 'figure' ? figureContent(block, text) : text;
    if (!preview) article.append(blockTools(block, index, content));
    article.append(content);
    paper.append(article);
  });
  if (!preview) paper.append(dropLine(current.blocks.length));
  markComments(paper, current.comments ?? []);
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
  try {
    const { pages: index, trash: deleted } = await request('/api/notebook/pages');
    const page = index.find((item) => item.id === hash) ?? index.find((item) => item.section === 'results') ?? index[0];
    if (page) await openPage(page.id);
    else {
      showEmptyWorkspace();
      renderNavigation(index, deleted);
    }
  } catch (error) {
    if (error.status === 401) showAccess();
    else {
      $('workspace').hidden = false;
      $('access').hidden = true;
      pageLoadFailed = true;
      updateControls();
      status(error.message, true);
    }
  }
  if (!$('workspace').hidden) await refreshLibrary();
}

for (const section of sections) section.querySelector('.new-page').addEventListener('click', () => showCreate(section));
$('create-page-form').addEventListener('submit', async (event) => {
  event.preventDefault();
  if (loading || organizing || current?.renaming || current?.postingComment) return;
  const id = await organizePages('create', '/api/notebook/pages', { section: createSection, title: $('new-page-title').value });
  if (id) await openPage(id);
});
$('create-page-cancel').addEventListener('click', () => {
  closePageForm('create');
  sections.find((section) => section.dataset.section === createSection).querySelector('.new-page').focus();
});
$('move-page').addEventListener('click', showMove);
$('move-section').addEventListener('change', movePositions);
$('move-page-form').addEventListener('submit', async (event) => {
  event.preventDefault();
  if (!current || loading || organizing || current.renaming) return;
  if (await organizePages('move', '/api/notebook/move', {
    page: current.id, section: $('move-section').value, after: $('move-position').value || null,
  })) $('move-page').focus();
});
$('move-page-cancel').addEventListener('click', () => {
  closePageForm('move');
  $('move-page').focus();
});
$('delete-page').addEventListener('click', () => void changeTrash(current.id, 'delete'));
$('rename').addEventListener('click', () => {
  closePageForm('create');
  closePageForm('move');
  $('rename-form').hidden = false;
  $('title-input').value = $('page-title').textContent;
  $('title-input').select();
  $('title-input').focus();
});
$('rename-form').addEventListener('submit', renamePage);
$('rename-cancel').addEventListener('click', () => {
  closeRename();
  $('rename').focus();
});
$('comments-jump').addEventListener('click', () => $('page-comments').scrollIntoView({ block: 'start' }));
$('comment-selection').addEventListener('pointerdown', (event) => event.preventDefault());
$('comment-selection').addEventListener('click', commentOnSelection);
$('comment-cancel').addEventListener('click', () => {
  commentDrafts.delete(current.id);
  showCommentDraft();
  commentStatus();
});
$('comments-refresh').addEventListener('click', () => void loadComments(current));
$('comment-form').addEventListener('submit', postComment);
$('comment-text').addEventListener('input', (event) => {
  const draft = commentDrafts.get(current.id);
  if (draft) draft.text = event.target.value;
});
document.addEventListener('selectionchange', () => {
  if (preview && current) offerSelection(previewAnchor(paper));
});
function setPreview(value) {
  if (preview === value) return;
  preview = value;
  $('preview').setAttribute('aria-pressed', String(preview));
  $('preview').textContent = preview ? 'Edit' : 'View';
  paper.classList.toggle('preview', preview);
  $('workspace').classList.toggle('viewing', preview);
  if (preview) {
    closeRename();
    closePageForm('create');
    closePageForm('move');
  }
  renderPage();
  updateControls();
}
$('preview').addEventListener('click', () => setPreview(!preview));
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
      await loadComments(current);
    } else await start();
  } catch (error) {
    $('login-error').textContent = error.message;
  }
});
window.addEventListener('beforeunload', (event) => {
  if (!current?.dirty && ![...commentDrafts.values()].some((draft) => draft.text.trim())) return;
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
