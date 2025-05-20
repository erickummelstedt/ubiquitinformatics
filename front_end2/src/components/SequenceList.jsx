import React from 'react';
import Visualizer from './Visualizer';

const MOL2_URL = '/his-GG-1ubq-1-(<K48_1ubq-2-()>).mol2';

// Utility to convert ^* to superscript and _* to subscript
const superscriptMap = {
  '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
  '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
  'a': 'ᵃ', 'b': 'ᵇ', 'c': 'ᶜ', 'd': 'ᵈ', 'e': 'ᵉ',
  'f': 'ᶠ', 'g': 'ᵍ', 'h': 'ʰ', 'i': 'ⁱ', 'j': 'ʲ',
  'k': 'ᵏ', 'l': 'ˡ', 'm': 'ᵐ', 'n': 'ⁿ', 'o': 'ᵒ',
  'p': 'ᵖ', 'r': 'ʳ', 's': 'ˢ', 't': 'ᵗ', 'u': 'ᵘ',
  'v': 'ᵛ', 'w': 'ʷ', 'x': 'ˣ', 'y': 'ʸ', 'z': 'ᶻ',
  'A': 'ᴬ', 'B': 'ᴮ', 'D': 'ᴰ', 'E': 'ᴱ', 'G': 'ᴳ',
  'H': 'ᴴ', 'I': 'ᴵ', 'J': 'ᴶ', 'K': 'ᴷ', 'L': 'ᴸ',
  'M': 'ᴹ', 'N': 'ᴺ', 'O': 'ᴼ', 'P': 'ᴾ', 'R': 'ᴿ',
  'T': 'ᵀ', 'U': 'ᵁ', 'V': 'ⱽ', 'W': 'ᵂ'
};
const subscriptMap = {
  '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄',
  '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'
};
function formatSuperSub(text) {
  // Replace ^* with superscript
  text = text.replace(/\^([A-Za-z0-9])/g, (m, c) => superscriptMap[c] || c);
  // Replace _* with subscript
  text = text.replace(/_([0-9])/g, (m, c) => subscriptMap[c] || c);
  // Replace superscript followed by space and number (optionally before \n or end) with bold number
  text = text.replace(/([\u00B2-\u00B3\u00B9\u2070-\u209F]) (\d+)(?=\n|$)/g, (m, sup, num) => `${sup} <span style='font-weight:900'>${num}</span>`);
  // Replace \n with <br> for new lines
  text = text.replace(/\n/g, '<br>');
  return text;
}

const SequenceList = ({ sequences, onSelect }) => {
  return (
    <div style={{ width: '100%' }}>
      <div style={{ maxHeight: '200px', overflowY: 'auto', marginBottom: '1rem', border: '1px solid #eee', borderRadius: '4px' }}>
        <ul style={{ listStyle: 'none', padding: 0, margin: 0 }}>
          {sequences.length === 0 && <li style={{ padding: '0.5rem' }}>No sequences found.</li>}
          {sequences.map((seqArr, idx) => {
            // Only show Acceptor and step keys
            const filtered = seqArr.filter(({ key }) => key.toLowerCase().includes('acceptor') || key.toLowerCase().includes('step'));
            return (
              <li
                key={idx}
                onClick={() => onSelect(seqArr)}
                style={{ cursor: 'pointer', padding: '0.5rem', borderBottom: '1px solid #f0f0f0' }}
              >
                <table style={{ borderCollapse: 'collapse', margin: '0 auto', width: 'auto', tableLayout: 'fixed' }}>
                  <thead>
                    <tr style={{ background: '#f7f7fa' }}>
                      {filtered.map(({ key }) => (
                        <th key={key} style={{
                          width: 140,
                          height: 48,
                          padding: '0 12px',
                          background: 'inherit',
                          fontWeight: 700,
                          fontSize: 14,
                          color: '#222',
                          textAlign: 'center',
                          whiteSpace: 'nowrap',
                          verticalAlign: 'middle',
                          border: 'none',
                        }} dangerouslySetInnerHTML={{ __html: formatSuperSub(key) }} />
                      ))}
                    </tr>
                    <tr>
                      <td colSpan={filtered.length} style={{ height: 0, padding: 0 }}>
                        <div style={{ borderBottom: '2px solid #bbb', width: '100%' }} />
                      </td>
                    </tr>
                  </thead>
                  <tbody>
                    <tr>
                      {filtered.map(({ key, value }) => (
                        <td key={key} style={{
                          width: 140,
                          height: 48,
                          padding: '0 12px',
                          background: '#fff',
                          fontWeight: 400, // set to normal so bold span is visible
                          fontSize: 15,
                          color: '#222',
                          textAlign: 'center',
                          whiteSpace: 'nowrap',
                          verticalAlign: 'middle',
                          border: 'none',
                        }} dangerouslySetInnerHTML={{ __html: formatSuperSub(String(value)) }} />
                      ))}
                    </tr>
                  </tbody>
                </table>
              </li>
            );
          })}
        </ul>
      </div>
      {/* Visualizer bar: 6 small 3Dmol.js viewers */}
      <div style={{ display: 'flex', flexDirection: 'row', justifyContent: 'center', alignItems: 'center', gap: 0, height: 180, width: '90vw', marginLeft: 'calc(-45vw + 50%)', marginTop: 16, background: '#181818', border: '1px solid #ccc', borderRadius: 50, padding: 16 }}>
        {[...Array(6)].map((_, i) => (
          <div key={i} style={{ height: 200, width: 200, position: 'relative', display: 'flex', alignItems: 'center', justifyContent: 'center', marginLeft: i !== 0 ? 16 : 0 }}>
            <div style={{ width: 200, height: 200, position: 'relative' }}>
              <Visualizer mol2Url={MOL2_URL} />
            </div>
          </div>
        ))}
      </div>
    </div>
  );
};

export default SequenceList;
