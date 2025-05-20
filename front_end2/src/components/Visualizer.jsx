import React, { useRef, useEffect } from 'react';

const Visualizer = ({ mol2Url }) => {
  const viewerRef = useRef(null);
  useEffect(() => {
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
      element.style.width = '100%';
      element.style.height = '100%';
      element.style.minHeight = '0';
      element.style.position = 'absolute';
      element.style.top = '0';
      element.style.left = '0';
      element.style.right = '0';
      element.style.bottom = '0';
      element.style.overflow = 'hidden';
      element.style.background = '#232323';
      element.style.border = '2px solid #333';
      const config = { backgroundColor: '#232323' };
      const viewer = window.$3Dmol.createViewer(element, config);
      window.fetch(mol2Url)
        .then(res => res.text())
        .then(mol2str => {
          viewer.addModel(mol2str, 'mol2');
          viewer.setStyle({}, { cartoon: { color: 'spectrum' }, stick: {} });
          viewer.zoomTo();
          viewer.render();
        });
    }
  }, [mol2Url]);
  return (
    <div style={{ position: 'relative', width: '100%', height: '100%', minHeight: 0, flex: 1, overflow: 'hidden', display: 'flex' }}>
      <div ref={viewerRef} style={{ width: '100%', height: '100%', minHeight: 0, position: 'absolute', top: 0, left: 0, right: 0, bottom: 0, overflow: 'hidden', background: '#232323', border: '2px solid #333' }} />
    </div>
  );
};

export default Visualizer;
