import React, { useState } from 'react';

const panel_width = 600;
const panel_height = 400;
const canvas_width = panel_width - 30;
const canvas_height = panel_height - 30;

const basePanelStyle = {
  backgroundColor: '#181818',
  padding: '10px',
  border: '1px solid',
  borderColor: '#444444',
  width: panel_width + 'px',
  height: panel_height + 'px',
  overflow: 'auto',
  boxSizing: 'border-box',
  color: '#E0E0E0',
  cursor: 'pointer',
  borderRadius: '16px',
  marginBottom: '1rem',
};

const SequenceFilter = ({ value, onChange, onPanelSelect }) => {
  // Panel component from main.js
  const Panel = ({ children }) => {
    const [hover, setHover] = useState(false);
    const style = {
      ...basePanelStyle,
      borderColor: hover ? '#888888' : '#444444',
    };
    return (
      <div
        style={style}
        onMouseEnter={() => setHover(true)}
        onMouseLeave={() => setHover(false)}
      >
        {children}
      </div>
    );
  };

  // GraphPanel placeholder (canvas logic can be added later)
  const GraphPanel = () => (
    <canvas
      width={canvas_width}
      height={canvas_height}
      style={{ display: 'block', margin: 'auto', backgroundColor: '#1c1c1c', border: '1px solid #888', width: canvas_width + 'px', height: canvas_height + 'px' }}
    />
  );

  return (
    <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center' }}>
      <input
        type="text"
        placeholder="Filter sequences..."
        value={value}
        onChange={onChange}
        style={{ width: '100%', maxWidth: panel_width + 'px', padding: '0.5rem', fontSize: '1rem', marginBottom: '1rem' }}
      />
      <Panel>
        <GraphPanel />
      </Panel>
      <div style={{ marginTop: '1rem', color: '#888', fontSize: '0.9rem', textAlign: 'center' }}>
        {/* Placeholder for connection to output selection and visualizer */}
        <em>Panel selection and output connection logic will be implemented here.</em>
      </div>
    </div>
  );
};

export default SequenceFilter;
