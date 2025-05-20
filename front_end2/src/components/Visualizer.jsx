import React, { useEffect, useRef } from 'react';

const MOL2_URL = process.env.NODE_ENV === 'development'
  ? '/his-GG-1ubq-1-(<K48_1ubq-2-()>).mol2'
  : 'his-GG-1ubq-1-(<K48_1ubq-2-()>).mol2';

const Visualizer = () => {
  const viewerRef = useRef(null);
  useEffect(() => {
    // Dynamically load 3Dmol.js if not already loaded
    if (!window.$3Dmol) {
      const script = document.createElement('script');
      script.src = 'https://3dmol.org/build/3Dmol-min.js';
      script.async = true;
      script.onload = () => initViewer();
      document.body.appendChild(script);
    } else {
      initViewer();
    }
    function initViewer() {
      if (!viewerRef.current) return;
      viewerRef.current.innerHTML = '';
      const element = viewerRef.current;
      const config = { backgroundColor: '#232323' };
      const viewer = window.$3Dmol.createViewer(element, config);
      window.fetch(MOL2_URL)
        .then(res => res.text())
        .then(mol2str => {
          viewer.addModel(mol2str, 'mol2');
          viewer.setStyle({}, { cartoon: { color: 'spectrum' }, stick: {} });
          viewer.zoomTo();
          viewer.render();
        });
    }
  }, []);
  return (
    <div style={{ width: '100%', height: '100%', minHeight: 0 }}>
      <div ref={viewerRef} style={{ width: '100%', height: '100%', minHeight: 0, flex: 1 }} />
    </div>
  );
};

export default Visualizer;
