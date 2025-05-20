import React, { useState, useRef, useEffect } from 'react';

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
  overflow: 'hidden', // prevent scrolling
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

  // --- Game-like Scaffold Panel ---
  const GameScaffoldPanel = () => {
    const canvasRef = useRef(null);
    const [refreshKey, setRefreshKey] = useState(0); // for refresh button

    // Scaffold state
    const RADIUS = 20;
    const LIGHT_GRAY = '#dddddd';
    const GRAY = '#aaaaaa';
    const BLACK = '#000000';
    // Panel dimensions
    const PANEL_WIDTH = canvas_width;
    const PANEL_HEIGHT = canvas_height;

    // Centered initial nodes for a 600x400 panel, shifted down by 50px
    const initialNodes = [
      { x: 300, y: 300, clicks: 0 },
      { x: 250, y: 250, clicks: 0 }, { x: 350, y: 250, clicks: 0 },
      { x: 200, y: 200, clicks: 0 }, { x: 300, y: 200, clicks: 0 }, { x: 400, y: 200, clicks: 0 },
      { x: 150, y: 150, clicks: 0 }, { x: 250, y: 150, clicks: 0 }, { x: 350, y: 150, clicks: 0 }, { x: 450, y: 150, clicks: 0 },
      { x: 100, y: 100, clicks: 0 }, { x: 200, y: 100, clicks: 0 }, { x: 300, y: 100, clicks: 0 }, { x: 400, y: 100, clicks: 0 }, { x: 500, y: 100, clicks: 0 }
    ];
    const edges = [
      [0, 1], [0, 2],
      [1, 3], [1, 4],
      [2, 4], [2, 5],
      [3, 6], [3, 7],
      [4, 7], [4, 8],
      [5, 8], [5, 9],
      [6, 10], [6, 11],
      [7, 11], [7, 12],
      [8, 12], [8, 13],
      [9, 13], [9, 14]
    ];

    // State for arrows, nodes, etc.
    const [nodes, setNodes] = useState(() => {
      const n = initialNodes.map(node => ({ ...node }));
      n[0].clicks = 1;
      return n;
    });
    const [arrows, setArrows] = useState([]);
    const [clickedNodes, setClickedNodes] = useState(new Set([0]));
    const [lastNodeIndex, setLastNodeIndex] = useState(0);
    const maxClicks = 5;
    const [refreshBtnRect, setRefreshBtnRect] = useState(null);

    // refs for nodes, arrows, and clickedNodes
    const nodesRef = useRef(nodes);
    const arrowsRef = useRef(arrows);
    const clickedNodesRef = useRef(clickedNodes);
    useEffect(() => { nodesRef.current = nodes; }, [nodes]);
    useEffect(() => { arrowsRef.current = arrows; }, [arrows]);
    useEffect(() => { clickedNodesRef.current = clickedNodes; }, [clickedNodes]);

    // Draw the scaffold
    useEffect(() => {
      const canvas = canvasRef.current;
      if (!canvas) return;
      const ctx = canvas.getContext('2d');
      // Set canvas to panel size
      canvas.width = PANEL_WIDTH;
      canvas.height = PANEL_HEIGHT;
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      // Remove BASE_WIDTH/BASE_HEIGHT scaling logic
      // Scaffold bounds (for box)
      const padding = 0;
      const xs = nodes.map(n => n.x);
      const ys = nodes.map(n => n.y);
      const minX = Math.min(...xs) - RADIUS - padding;
      const maxX = Math.max(...xs) + RADIUS + padding;
      const minY = Math.min(...ys) - RADIUS - padding;
      const maxY = Math.max(...ys) + RADIUS + padding;
      // Draw background box (now fills the panel)
      const boxX = 0;
      const boxY = 0;
      const boxW = PANEL_WIDTH;
      const boxH = PANEL_HEIGHT;
      const radius = 20;
      ctx.fillStyle = '#101010';
      ctx.beginPath();
      ctx.moveTo(boxX + radius, boxY);
      ctx.arcTo(boxX + boxW, boxY, boxX + boxW, boxY + boxH, radius);
      ctx.arcTo(boxX + boxW, boxY + boxH, boxX, boxY + boxH, radius);
      ctx.arcTo(boxX, boxY + boxH, boxX, boxY, radius);
      ctx.arcTo(boxX, boxY, boxX + boxW, boxY, radius);
      ctx.closePath();
      ctx.fill();
      // Draw refresh button (bottom right)
      const refreshBtn = {
        x: boxX + boxW - 100,
        y: boxY + boxH - 50,
        w: 80,
        h: 30,
        r: 10
      };
      ctx.fillStyle = '#222';
      ctx.beginPath();
      ctx.moveTo(refreshBtn.x + refreshBtn.r, refreshBtn.y);
      ctx.arcTo(refreshBtn.x + refreshBtn.w, refreshBtn.y, refreshBtn.x + refreshBtn.w, refreshBtn.y + refreshBtn.h, refreshBtn.r);
      ctx.arcTo(refreshBtn.x + refreshBtn.w, refreshBtn.y + refreshBtn.h, refreshBtn.x, refreshBtn.y + refreshBtn.h, refreshBtn.r);
      ctx.arcTo(refreshBtn.x, refreshBtn.y + refreshBtn.h, refreshBtn.x, refreshBtn.y, refreshBtn.r);
      ctx.arcTo(refreshBtn.x, refreshBtn.y, refreshBtn.x + refreshBtn.w, refreshBtn.y, refreshBtn.r);
      ctx.closePath();
      ctx.fill();
      ctx.strokeStyle = '#555';
      ctx.lineWidth = 1.5;
      ctx.stroke();
      ctx.fillStyle = '#fff';
      ctx.font = `16px 'Helvetica Neue', sans-serif`;
      ctx.fillText('Refresh', refreshBtn.x + 10, refreshBtn.y + 20);
      setRefreshBtnRect(refreshBtn);
      // Outer/inner border
      ctx.strokeStyle = '#444';
      ctx.lineWidth = 4;
      ctx.stroke();
      ctx.strokeStyle = '#666';
      ctx.lineWidth = 2;
      ctx.stroke();
      // Draw edges (no scaling/offset)
      ctx.strokeStyle = GRAY;
      ctx.lineWidth = 2;
      edges.forEach(([i1, i2]) => {
        const { x: x1, y: y1 } = nodes[i1];
        const { x: x2, y: y2 } = nodes[i2];
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
      });
      // Draw arrows (no scaling/offset)
      arrows.forEach(({ from, to, color }) => {
        const { x: x1, y: y1 } = nodes[from];
        const { x: x2, y: y2 } = nodes[to];
        const dx = x2 - x1;
        const dy = y2 - y1;
        const length = Math.sqrt(dx * dx + dy * dy);
        const offset = 5;
        const offsetXperp = -dy / length * offset;
        const offsetYperp = dx / length * offset;
        const direction = color === 'blue' ? 1 : -1;
        const sx1 = x1 + direction * offsetXperp;
        const sy1 = y1 + direction * offsetYperp;
        const sx2 = x2 + direction * offsetXperp;
        const sy2 = y2 + direction * offsetYperp;
        ctx.strokeStyle = color;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(sx1, sy1);
        ctx.lineTo(sx2, sy2);
        ctx.stroke();
        // Arrowhead
        const midX = (sx1 + sx2) / 2;
        const midY = (sy1 + sy2) / 2;
        const angle = Math.atan2(sy2 - sy1, sx2 - sx1);
        const arrowLength = 10;
        ctx.fillStyle = color;
        ctx.beginPath();
        ctx.moveTo(midX, midY);
        ctx.lineTo(
          midX - arrowLength * Math.cos(angle - Math.PI / 6),
          midY - arrowLength * Math.sin(angle - Math.PI / 6)
        );
        ctx.lineTo(
          midX - arrowLength * Math.cos(angle + Math.PI / 6),
          midY - arrowLength * Math.sin(angle + Math.PI / 6)
        );
        ctx.closePath();
        ctx.fill();
      });
      // Draw nodes (no scaling/offset)
      nodes.forEach(({ x, y, clicks }, i) => {
        if (clicks === 0) ctx.fillStyle = LIGHT_GRAY;
        else if (clicks === 1) ctx.fillStyle = '#ff9999';
        else ctx.fillStyle = '#cc6666';
        ctx.beginPath();
        ctx.arc(x, y, RADIUS, 0, Math.PI * 2);
        ctx.fill();
        ctx.strokeStyle = BLACK;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(x, y, RADIUS, 0, Math.PI * 2);
        ctx.stroke();
      });
    }, [nodes, arrows, refreshKey]);

    // --- FULL INTERACTIVE CLICK HANDLER WITH ARROWS (main copy.js logic) ---
    useEffect(() => {
      const canvas = canvasRef.current;
      if (!canvas) return;
      const handleMouseDown = (e) => {
        const rect = canvas.getBoundingClientRect();
        const x = (e.clientX - rect.left);
        const y = (e.clientY - rect.top);
        // Check refresh button
        if (refreshBtnRect) {
          const { x: bx, y: by, w: bw, h: bh } = refreshBtnRect;
          if (x >= bx && x <= bx + bw && y >= by && y <= by + bh) {
            setNodes(() => {
              const n = initialNodes.map(node => ({ ...node }));
              n[0].clicks = 1;
              return n;
            });
            setArrows([]);
            setClickedNodes(new Set([0]));
            setLastNodeIndex(0);
            setRefreshKey(k => k + 1);
            return;
          }
        }
        for (let i = 0; i < nodesRef.current.length; i++) {
          const node = nodesRef.current[i];
          const drawX = node.x;
          const drawY = node.y;
          const dx = drawX - x;
          const dy = drawY - y;
          if (Math.sqrt(dx * dx + dy * dy) <= RADIUS) {
            let isAllowed = false;
            let proposedPrev = null;
            Array.from(clickedNodesRef.current).forEach(prevIndex => {
              const isAdjacent = edges.some(([a, b]) =>
                (a === prevIndex && b === i) || (b === prevIndex && a === i)
              );
              const isAbove = nodesRef.current[i].y < nodesRef.current[prevIndex].y;
              const lineExists = arrowsRef.current.some(({ from, to }) =>
                (from === prevIndex && to === i)
              );
              if (isAdjacent && isAbove && !lineExists) {
                isAllowed = true;
                proposedPrev = prevIndex;
              }
            });
            if (!isAllowed) continue;
            const totalClicks = nodesRef.current.reduce((sum, node) => sum + node.clicks, 0);
            if (totalClicks >= maxClicks) return;
            setNodes(prevNodes => prevNodes.map((n, idx) => idx === i ? { ...n, clicks: (n.clicks || 0) + 1 } : n));
            setClickedNodes(prev => new Set(prev).add(i));
            const prevIdx = proposedPrev !== null ? proposedPrev : i;
            if (prevIdx !== null && prevIdx !== i) {
              const prevNode = nodesRef.current[prevIdx];
              const newNode = nodesRef.current[i];
              const color = newNode.x < prevNode.x ? 'red' : 'blue';
              setArrows(prev => [...prev, { from: prevIdx, to: i, color }]);
            }
            setLastNodeIndex(i);
            return;
          }
        }
      };
      canvas.addEventListener('mousedown', handleMouseDown);
      return () => {
        canvas.removeEventListener('mousedown', handleMouseDown);
      };
    }, [nodes, arrows, clickedNodes]);

    return (
      <div style={{ position: 'relative', width: '100%', height: '100%' }}>
        <canvas ref={canvasRef} style={{ width: '100%', height: '100%' }} />
      </div>
    );
  };

  // --- Visualizer (3Dmol.js) Panel ---
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
        // Strictly enforce size and containment
        element.style.width = '100%';
        element.style.height = '100%';
        element.style.minHeight = '0';
        element.style.position = 'absolute'; // absolute for full fill
        element.style.top = '0';
        element.style.left = '0';
        element.style.right = '0';
        element.style.bottom = '0';
        element.style.overflow = 'hidden';
        element.style.background = '#232323';
        element.style.border = '2px solid #333';
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
    // Use a relative parent and absolute child for strict containment
    return (
      <div style={{ position: 'relative', width: '100%', height: '100%', minHeight: 0, flex: 1, overflow: 'hidden', display: 'flex' }}>
        <div ref={viewerRef} style={{ width: '100%', height: '100%', minHeight: 0, position: 'absolute', top: 0, left: 0, right: 0, bottom: 0, overflow: 'hidden', background: '#232323', border: '2px solid #333' }} />
      </div>
    );
  };

  return (
    <div style={{ display: 'flex', flexDirection: 'row', alignItems: 'flex-start', justifyContent: 'center', width: '100%' }}>
      {/* Left/Main Panel */}
      <Panel>
        <GameScaffoldPanel />
      </Panel>
      {/* Right Hand Panel with Visualizer */}
      <Panel>
        <div style={{ display: 'flex', flexDirection: 'column', height: '100%', width: '100%', position: 'relative', overflow: 'hidden' }}>
          <div style={{ flex: 1, width: '100%', minHeight: 0, position: 'relative', overflow: 'hidden' }}>
            <Visualizer />
          </div>
        </div>
      </Panel>
    </div>
  );
};

export default SequenceFilter;
