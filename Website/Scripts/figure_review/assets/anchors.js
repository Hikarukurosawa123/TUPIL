import { renderMarkdown } from './markdown.js';

export function blockId() {
  const bytes = crypto.getRandomValues(new Uint8Array(16));
  return `block-${Array.from(bytes, (byte) => byte.toString(16).padStart(2, '0')).join('')}`;
}

export function sourceAnchor(field) {
  const quote = field.value.slice(field.selectionStart, field.selectionEnd);
  return quote.trim() ? { block: field.dataset.block, mode: 'source', offset: field.selectionStart, quote } : null;
}

export function previewAnchor(paper) {
  const selection = document.getSelection();
  if (!selection?.rangeCount || selection.isCollapsed) return null;
  const range = selection.getRangeAt(0);
  const startNode = range.startContainer;
  const anchor = (startNode instanceof Element ? startNode : startNode.parentElement)?.closest('.markdown[data-block]');
  if (!anchor || !paper.contains(anchor) || !anchor.contains(range.endContainer)) return null;
  const prefix = new Range();
  prefix.selectNodeContents(anchor);
  prefix.setEnd(range.startContainer, range.startOffset);
  const quote = range.toString();
  return quote.trim() ? { block: anchor.dataset.block, mode: 'preview', offset: prefix.toString().length, quote } : null;
}

function locate(text, quote, previous) {
  if (!quote) return null;
  if (previous !== null && text.slice(previous, previous + quote.length) === quote) return [previous, previous + quote.length];
  const start = text.indexOf(quote);
  return start >= 0 && text.indexOf(quote, start + 1) < 0 ? [start, start + quote.length] : null;
}

function toRange(element, [start, end]) {
  const walker = document.createTreeWalker(element, NodeFilter.SHOW_TEXT);
  const range = new Range();
  let offset = 0;
  let started = false;
  for (let text = walker.nextNode(); text; text = walker.nextNode()) {
    const next = offset + text.length;
    if (!started && start < next) {
      range.setStart(text, start - offset);
      started = true;
    }
    if (started && end <= next) {
      range.setEnd(text, end - offset);
      return range;
    }
    offset = next;
  }
  return null;
}

function previewRange(element, anchor) {
  const text = element.textContent;
  const isPreview = anchor.mode === 'preview';
  const offsets = locate(text, anchor.quote, isPreview ? anchor.offset : null)
    ?? (!isPreview ? locate(text, renderMarkdown(anchor.quote).textContent, null) : null);
  return offsets && toRange(element, offsets);
}

export function markComments(paper, comments) {
  if (!globalThis.CSS?.highlights || !globalThis.Highlight) return;
  CSS.highlights.delete('comment-selection');
  const highlight = new Highlight();
  const blocks = new Map([...paper.querySelectorAll('.markdown[data-block]')].map((item) => [item.dataset.block, item]));
  for (const { anchor } of comments) {
    const block = anchor && blocks.get(anchor.block);
    const range = block && previewRange(block, anchor);
    if (range) highlight.add(range);
  }
  if (highlight.size) CSS.highlights.set('comment-selection', highlight);
}

export function revealAnchor(paper, anchor) {
  const target = paper.querySelector(`[data-block="${CSS.escape(anchor.block)}"]`);
  if (!target) return false;
  if (anchor.mode === 'source') {
    const offsets = locate(target.value, anchor.quote, anchor.offset);
    if (!offsets) return false;
    target.scrollIntoView({ block: 'center' });
    target.focus({ preventScroll: true });
    target.setSelectionRange(...offsets);
  } else {
    const range = previewRange(target, anchor);
    if (!range) return false;
    const selection = document.getSelection();
    selection.removeAllRanges();
    selection.addRange(range);
    const start = range.startContainer.parentElement;
    start.scrollIntoView({ block: 'center' });
  }
  return true;
}
