import React from 'react';
import Visualizer from './Visualizer';

const MOL2_URL = '/his-GG-1ubq-1-(<K48_1ubq-2-()>).mol2';

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
                    <tr>
                      {filtered.map(({ key }) => (
                        <th key={key} style={{
                          minWidth: 80,
                          maxWidth: 160,
                          padding: '6px 12px',
                          background: '#f7f7fa',
                          border: '1px solid #d0d0e0',
                          fontWeight: 700,
                          fontSize: 13,
                          color: '#222',
                          textAlign: 'center',
                          whiteSpace: 'nowrap',
                        }}>{key}</th>
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    <tr>
                      {filtered.map(({ key, value }) => (
                        <td key={key} style={{
                          minWidth: 80,
                          maxWidth: 160,
                          padding: '6px 12px',
                          background: '#fff',
                          border: '1px solid #d0d0e0',
                          fontWeight: 500,
                          fontSize: 14,
                          color: '#222',
                          textAlign: 'center',
                          whiteSpace: 'nowrap',
                        }}>{String(value)}</td>
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
