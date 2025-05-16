import React from 'react';

const SequenceList = ({ sequences, onSelect }) => {
  return (
    <div style={{ maxHeight: '200px', overflowY: 'auto', marginBottom: '1rem', border: '1px solid #eee', borderRadius: '4px' }}>
      <ul style={{ listStyle: 'none', padding: 0, margin: 0 }}>
        {sequences.length === 0 && <li style={{ padding: '0.5rem' }}>No sequences found.</li>}
        {sequences.map((seq, idx) => (
          <li
            key={idx}
            onClick={() => onSelect(seq)}
            style={{ cursor: 'pointer', padding: '0.5rem', borderBottom: '1px solid #f0f0f0' }}
          >
            {Object.keys(seq).join(', ')}
          </li>
        ))}
      </ul>
    </div>
  );
};

export default SequenceList;
