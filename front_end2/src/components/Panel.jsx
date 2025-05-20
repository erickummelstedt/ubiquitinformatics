import React from 'react';

const Panel = ({ children, style = {}, ...props }) => {
  return (
    <div
      style={{
        backgroundColor: '#181818',
        padding: '10px',
        border: '1px solid #444444',
        width: 600,
        height: 400,
        overflow: 'hidden',
        boxSizing: 'border-box',
        color: '#E0E0E0',
        cursor: 'pointer',
        borderRadius: '16px',
        marginBottom: '1rem',
        ...style,
      }}
      {...props}
    >
      {children}
    </div>
  );
};

export default Panel;
