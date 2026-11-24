const tokens = /(`[^`\n]+`|\*\*[^*\n]+\*\*|__[^_\n]+__|\*[^*\n]+\*|\b_[^_\n]+_\b|\[[^\]\n]+\]\([^\s)]+\))/g;
const listPattern = /^\s*(\d+\.|[-+*])\s+(.*)$/;
const headingPattern = /^(#{1,6})\s+(.*)$/;
const quotePattern = /^>\s?(.*)$/;
const fencePattern = /^\s*```[\w-]*\s*$/;
const rulePattern = /^\s*(?:-{3,}|\*{3,}|_{3,})\s*$/;

function inline(element, text) {
  for (const [index, part] of text.split(tokens).entries()) {
    if (index % 2 === 0) {
      element.append(part);
      continue;
    }
    let child;
    if (part.startsWith('[')) {
      const split = part.lastIndexOf('](');
      const target = part.slice(split + 2, -1);
      if (!/^(https?:\/\/|mailto:)/i.test(target)) {
        element.append(part);
        continue;
      }
      child = document.createElement('a');
      child.href = target;
      child.target = '_blank';
      child.rel = 'noopener noreferrer';
      child.textContent = part.slice(1, split);
    } else {
      const markerLength = part.startsWith('**') || part.startsWith('__') ? 2 : 1;
      child = document.createElement(part[0] === '`' ? 'code' : markerLength === 2 ? 'strong' : 'em');
      child.textContent = part.slice(markerLength, -markerLength);
    }
    element.append(child);
  }
}

export function renderMarkdown(text) {
  const result = document.createElement('div');
  result.className = 'markdown';
  const lines = text.replace(/\r\n?/g, '\n').split('\n');
  let paragraph;
  let list;
  let quote;
  let code;
  for (const line of lines) {
    if (code) {
      if (/^\s*```\s*$/.test(line)) code = null;
      else code.textContent += `${line}\n`;
      continue;
    }
    const heading = line.match(headingPattern);
    const item = line.match(listPattern);
    const quoted = line.match(quotePattern);
    if (!item) list = null;
    if (!quoted) quote = null;
    if (fencePattern.test(line)) {
      paragraph = null;
      const pre = document.createElement('pre');
      code = document.createElement('code');
      pre.append(code);
      result.append(pre);
    } else if (heading || item || quoted || rulePattern.test(line) || !line.trim()) {
      paragraph = null;
      if (heading) {
        const title = document.createElement(`h${heading[1].length}`);
        inline(title, heading[2]);
        result.append(title);
      } else if (rulePattern.test(line)) {
        result.append(document.createElement('hr'));
      } else if (item) {
        const tag = item[1].endsWith('.') ? 'ol' : 'ul';
        if (list?.localName !== tag) {
          list = document.createElement(tag);
          if (tag === 'ol') list.start = Number.parseInt(item[1], 10);
          result.append(list);
        }
        const entry = document.createElement('li');
        inline(entry, item[2]);
        list.append(entry);
      } else if (quoted) {
        if (!quote) {
          quote = document.createElement('blockquote');
          result.append(quote);
        } else quote.append(document.createElement('br'));
        inline(quote, quoted[1]);
      }
    } else {
      if (!paragraph) {
        paragraph = document.createElement('p');
        result.append(paragraph);
      } else paragraph.append(document.createElement('br'));
      inline(paragraph, line);
    }
  }
  return result;
}
