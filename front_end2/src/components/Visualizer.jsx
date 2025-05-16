import React from 'react';

const Visualizer = ({ data }) => {
  if (!data) {
    return <div style={{ minHeight: '120px', background: '#f5f5f5', borderRadius: '4px', padding: '1rem', gridColumn: '1 / -1' }}>No data selected</div>;
  }

  // Display the first key/value pairs in a readable way
  return (
    <div style={{ minHeight: '120px', background: '#f5f5f5', borderRadius: '4px', padding: '1rem', gridColumn: '1 / -1' }}>
      <h2>Selected Reaction Sequence</h2>
      <table style={{ width: '100%', borderCollapse: 'collapse' }}>
        <tbody>
          {Object.entries(data).map(([key, value]) => (
            <tr key={key}>
              <td style={{ fontWeight: 'bold', padding: '0.25rem 0.5rem', borderBottom: '1px solid #eee' }}>{key}</td>
              <td style={{ padding: '0.25rem 0.5rem', borderBottom: '1px solid #eee', wordBreak: 'break-all' }}>{value}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
};

export default Visualizer;
