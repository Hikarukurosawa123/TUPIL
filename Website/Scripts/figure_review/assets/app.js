const byId = (id) => document.getElementById(id);
const elements = Object.fromEntries([
  'access', 'login', 'code', 'login-error', 'workspace', 'refresh', 'search', 'figures',
  'empty', 'empty-title', 'empty-message', 'viewer', 'folder', 'figure-title', 'open-figure',
  'surface', 'media', 'pins', 'status', 'add-note', 'note-help', 'note-form', 'note-text', 'cancel', 'note-list',
].map((id) => [id, byId(id)]));

const state = { figures: [], figure: null, notes: [], selected: null, draft: null, drag: null, pending: false };
const nameParts = (key) => key.split('/');
const figureQuery = (figure) => new URLSearchParams({ figure });

async function request(path, body) {
  const options = body ? {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body),
  } : {};
  const response = await fetch(path, options);
  const result = await response.json();
  if (!response.ok) {
    const error = new Error(result.error);
    error.status = response.status;
    throw error;
  }
  return result;
}

function status(message = '', error = false) {
  elements.status.textContent = message;
  elements.status.classList.toggle('error', error);
}

function empty(title, message) {
  elements.empty.hidden = false;
  elements.viewer.hidden = true;
  elements['empty-title'].textContent = title;
  elements['empty-message'].textContent = message;
}

function renderLibrary() {
  const filter = elements.search.value.trim().toLowerCase();
  const matches = state.figures.filter((key) => key.toLowerCase().includes(filter));
  elements.figures.replaceChildren();
  for (const key of matches) {
    const parts = nameParts(key);
    const row = document.createElement('button');
    row.className = 'figure-row';
    row.textContent = parts.pop();
    row.title = key;
    if (parts.length) {
      const folder = document.createElement('span');
      folder.className = 'muted';
      folder.textContent = parts.join(' / ');
      row.append(folder);
    }
    if (key === state.figure) row.setAttribute('aria-current', 'true');
    row.addEventListener('click', () => selectFigure(key));
    elements.figures.append(row);
  }
  if (filter && !matches.length) {
    const message = document.createElement('p');
    message.className = 'muted';
    message.textContent = 'No matching figures.';
    elements.figures.append(message);
  }
}

function showAccess() {
  elements.access.hidden = false;
  elements.workspace.hidden = true;
  elements.refresh.hidden = true;
  elements.code.focus();
}

async function refresh() {
  if (state.pending || state.drag) return;
  elements.refresh.disabled = true;
  try {
    const { figures } = await request('/api/figures');
    elements.access.hidden = true;
    elements.workspace.hidden = false;
    elements.refresh.hidden = false;
    state.figures = figures;
    if (figures.includes(state.figure)) {
      renderLibrary();
      await loadNotes(state.figure);
      status('Figures and notes refreshed.');
    } else {
      await selectFigure(figures[0] ?? null);
    }
  } catch (error) {
    if (error.status === 401) {
      showAccess();
      return;
    }
    elements.workspace.hidden = false;
    elements.refresh.hidden = false;
    if (!state.figure) empty('Could not load figures', error.message);
    else status(error.message, true);
  } finally {
    elements.refresh.disabled = false;
  }
}

async function selectFigure(key) {
  if (state.pending || state.drag) return;
  if (state.draft?.mode !== 'move' && elements['note-text'].value.trim()
      && !window.confirm('Discard the unsaved note?')) return;
  state.figure = key;
  state.selected = null;
  state.notes = [];
  clearDraft();
  renderLibrary();
  status();
  if (!key) {
    elements.media.replaceChildren();
    empty('No figures yet', 'Add a figure to Google Drive’s Results_Main folder, then refresh.');
    return;
  }
  const parts = nameParts(key);
  const name = parts.pop();
  const url = `/figure?${figureQuery(key)}`;
  elements['figure-title'].textContent = name;
  elements.folder.textContent = parts.join(' / ');
  elements['open-figure'].href = url;
  const isPdf = name.toLowerCase().endsWith('.pdf');
  const preview = document.createElement(isPdf ? 'iframe' : 'img');
  if (isPdf) preview.title = name;
  else preview.alt = name;
  preview.addEventListener('error', () => {
    if (state.figure === key) status('The figure could not be opened. Check Google Drive and refresh.', true);
  });
  preview.src = url;
  elements.media.replaceChildren(preview);
  elements.empty.hidden = true;
  elements.viewer.hidden = false;
  try {
    await loadNotes(key);
  } catch (error) {
    if (state.figure === key) status(error.message, true);
  }
}

async function loadNotes(key) {
  const { notes } = await request(`/api/notes?${figureQuery(key)}`);
  if (state.figure !== key || state.drag) return;
  state.notes = notes;
  renderNotes();
}

function position(pin, point) {
  pin.style.left = `${point.x * 100}%`;
  pin.style.top = `${point.y * 100}%`;
}

function pointAt(event) {
  const rect = elements.surface.getBoundingClientRect();
  const clamp = (value) => Math.max(0, Math.min(1, value));
  return { x: clamp((event.clientX - rect.left) / rect.width), y: clamp((event.clientY - rect.top) / rect.height) };
}

function button(label, action) {
  const control = document.createElement('button');
  control.type = 'button';
  control.className = 'button';
  control.textContent = label;
  control.addEventListener('click', action);
  return control;
}

function showNote(id) {
  state.selected = id;
  renderNotes();
  byId(`note-${id}`)?.scrollIntoView({ block: 'nearest' });
}

function makePin(note, index) {
  const pin = document.createElement('button');
  pin.type = 'button';
  pin.className = 'pin';
  pin.textContent = index + 1;
  pin.title = note.text;
  pin.setAttribute('aria-label', `Show note ${index + 1}`);
  pin.setAttribute('aria-pressed', String(state.selected === note.id));
  position(pin, note);
  pin.addEventListener('click', () => {
    if (!state.pending && !elements.surface.classList.contains('placing')) showNote(note.id);
  });
  pin.addEventListener('pointerdown', (event) => {
    if (event.button !== 0 || state.pending || state.draft || state.drag) return;
    pin.setPointerCapture(event.pointerId);
    state.drag = { pointerId: event.pointerId, startX: event.clientX, startY: event.clientY, moved: false };
  });
  pin.addEventListener('pointermove', (event) => {
    const drag = state.drag;
    if (drag?.pointerId !== event.pointerId) return;
    if (Math.hypot(event.clientX - drag.startX, event.clientY - drag.startY) > 4) drag.moved = true;
    if (drag.moved) position(pin, pointAt(event));
  });
  pin.addEventListener('pointerup', (event) => {
    const drag = state.drag;
    if (drag?.pointerId !== event.pointerId) return;
    state.drag = null;
    if (drag.moved) {
      position(pin, note);
      void save({ action: 'move', id: note.id, ...pointAt(event) }, 'Note moved.');
    }
  });
  pin.addEventListener('lostpointercapture', (event) => {
    if (state.drag?.pointerId === event.pointerId) {
      state.drag = null;
      position(pin, note);
    }
  });
  return pin;
}

function renderNotes() {
  const draft = state.draft;
  const placing = draft && !draft.point && draft.mode !== 'edit';
  elements.surface.classList.toggle('placing', Boolean(placing));
  elements.pins.replaceChildren();
  elements['note-list'].replaceChildren();
  elements['add-note'].disabled = !state.figure || Boolean(draft) || state.pending;
  elements['note-form'].hidden = !draft || draft.mode === 'move';
  elements['note-form'].querySelector('button[type="submit"], button:not([type])').disabled = Boolean(placing) || state.pending;
  let help = 'Select a figure to add notes.';
  if (state.figure) help = 'Click a pin to read its note. Drag a pin to move it.';
  if (draft?.mode === 'add') help = placing ? 'Click on or around the figure to place the note.' : 'Write your note and save it.';
  if (draft?.mode === 'edit') help = 'Edit the selected note.';
  if (draft?.mode === 'move') help = 'Click on or around the figure to move the pin. Press Escape to cancel.';
  elements['note-help'].textContent = help;
  state.notes.forEach((note, index) => {
    elements.pins.append(makePin(note, index));
    const item = document.createElement('li');
    item.id = `note-${note.id}`;
    item.className = `note-card${state.selected === note.id ? ' selected' : ''}`;
    const label = document.createElement('strong');
    label.textContent = `Note ${index + 1}`;
    const text = document.createElement('p');
    text.textContent = note.text;
    const actions = document.createElement('div');
    actions.className = 'actions';
    actions.append(
      button('Edit', () => startDraft('edit', note)),
      button('Move', () => startDraft('move', note)),
      button('Delete', () => {
        if (!state.pending && !state.draft && window.confirm('Delete this note?')) {
          void save({ action: 'delete', id: note.id }, 'Note deleted.');
        }
      }),
    );
    for (const control of actions.children) control.disabled = state.pending || Boolean(draft);
    item.append(label, text, actions);
    elements['note-list'].append(item);
  });
  if (draft?.mode === 'add' && draft.point) {
    const pin = document.createElement('span');
    pin.className = 'pin draft';
    pin.textContent = '+';
    pin.setAttribute('aria-hidden', 'true');
    position(pin, draft.point);
    elements.pins.append(pin);
  }
  if (draft?.mode === 'move') {
    const cancel = document.createElement('li');
    cancel.append(button('Cancel move', clearDraft));
    elements['note-list'].prepend(cancel);
  }
}

function startDraft(mode, note) {
  if (state.draft || state.pending) return;
  state.draft = { mode, id: note?.id, point: null };
  state.selected = note?.id ?? null;
  elements['note-text'].value = mode === 'edit' ? note.text : '';
  renderNotes();
  if (mode === 'edit') elements['note-text'].focus();
  else {
    elements.surface.tabIndex = 0;
    elements.surface.setAttribute('aria-label', 'Choose a position. Press Enter for the center, or click a location.');
    elements.surface.focus();
  }
}

function finishPlacement(point) {
  const draft = state.draft;
  if (!draft || draft.point || draft.mode === 'edit' || state.pending) return;
  if (draft.mode === 'move') {
    void save({ action: 'move', id: draft.id, ...point }, 'Note moved.');
  } else {
    draft.point = point;
    renderNotes();
    elements['note-text'].focus();
  }
}

function clearDraft() {
  state.draft = null;
  elements['note-text'].value = '';
  elements.surface.removeAttribute('tabindex');
  elements.surface.removeAttribute('aria-label');
  renderNotes();
}

async function save(payload, message) {
  if (state.pending) return;
  const figure = state.figure;
  state.pending = true;
  renderNotes();
  try {
    const { id } = await request('/api/notes', { ...payload, figure });
    state.selected = payload.action === 'delete' ? null : id;
    clearDraft();
    await loadNotes(figure);
    status(message);
  } catch (error) {
    if (error.status === 401) showAccess();
    else status(error.message, true);
  } finally {
    state.pending = false;
    renderNotes();
  }
}

elements.search.addEventListener('input', renderLibrary);
elements.refresh.addEventListener('click', refresh);
elements['add-note'].addEventListener('click', () => startDraft('add'));
elements.cancel.addEventListener('click', clearDraft);
elements.pins.addEventListener('click', (event) => finishPlacement(pointAt(event)));
elements.surface.addEventListener('keydown', (event) => {
  if (event.target === elements.surface && event.key === 'Enter') finishPlacement({ x: .5, y: .5 });
});
document.addEventListener('keydown', (event) => {
  if (event.key === 'Escape' && state.draft && !state.pending) clearDraft();
});
elements['note-form'].addEventListener('submit', (event) => {
  event.preventDefault();
  const draft = state.draft;
  if (!draft || (draft.mode === 'add' && !draft.point)) return;
  void save({ action: draft.mode, id: draft.id, text: elements['note-text'].value, ...draft.point }, 'Note saved.');
});
elements.login.addEventListener('submit', async (event) => {
  event.preventDefault();
  elements['login-error'].textContent = '';
  try {
    await request('/api/login', { code: elements.code.value });
    elements.login.reset();
    await refresh();
  } catch (error) {
    elements['login-error'].textContent = error.message;
  }
});

void refresh();
