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
    // Calculate scaling and centering
    const BASE_WIDTH = 900;
    const BASE_HEIGHT = 600;

    // Initial nodes and edges (static for now)
    const initialNodes = [
      { x: 400, y: 300, clicks: 0 },
      { x: 350, y: 250, clicks: 0 }, { x: 450, y: 250, clicks: 0 },
      { x: 300, y: 200, clicks: 0 }, { x: 400, y: 200, clicks: 0 }, { x: 500, y: 200, clicks: 0 },
      { x: 250, y: 150, clicks: 0 }, { x: 350, y: 150, clicks: 0 }, { x: 450, y: 150, clicks: 0 }, { x: 550, y: 150, clicks: 0 },
      { x: 200, y: 100, clicks: 0 }, { x: 300, y: 100, clicks: 0 }, { x: 400, y: 100, clicks: 0 }, { x: 500, y: 100, clicks: 0 }, { x: 600, y: 100, clicks: 0 }
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
      // Calculate scale to fit BASE_WIDTH x BASE_HEIGHT into panel
      const scale = Math.min(PANEL_WIDTH / BASE_WIDTH, PANEL_HEIGHT / BASE_HEIGHT);
      // Centering offsets
      const offsetX = (PANEL_WIDTH - BASE_WIDTH * scale) / 2;
      const offsetY = (PANEL_HEIGHT - BASE_HEIGHT * scale) / 2;

      // Scaffold bounds (for box)
      const padding = 30;
      const xs = nodes.map(n => n.x);
      const ys = nodes.map(n => n.y);
      const minX = Math.min(...xs) - RADIUS - padding;
      const maxX = Math.max(...xs) + RADIUS + padding;
      const minY = Math.min(...ys) - RADIUS - padding;
      const maxY = Math.max(...ys) + RADIUS + padding;

      // Draw background box (scaled and centered)
      const boxX = (minX * scale) + offsetX;
      const boxY = (minY * scale) + offsetY;
      const boxW = (maxX - minX) * scale;
      const boxH = (maxY - minY) * scale;
      const radius = 20 * scale;
      ctx.fillStyle = '#101010';
      ctx.beginPath();
      ctx.moveTo(boxX + radius, boxY);
      ctx.arcTo(boxX + boxW, boxY, boxX + boxW, boxY + boxH, radius);
      ctx.arcTo(boxX + boxW, boxY + boxH, boxX, boxY + boxH, radius);
      ctx.arcTo(boxX, boxY + boxH, boxX, boxY, radius);
      ctx.arcTo(boxX, boxY, boxX + boxW, boxY, radius);
      ctx.closePath();
      ctx.fill();

      // Draw refresh button (scaled and centered)
      const refreshBtn = {
        x: boxX + boxW - 100 * scale,
        y: boxY + boxH - 50 * scale,
        w: 80 * scale,
        h: 30 * scale,
        r: 10 * scale
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
      ctx.lineWidth = 1.5 * scale;
      ctx.stroke();
      ctx.fillStyle = '#fff';
      ctx.font = `${16 * scale}px 'Helvetica Neue', sans-serif`;
      ctx.fillText('Refresh', refreshBtn.x + 10 * scale, refreshBtn.y + 20 * scale);
      setRefreshBtnRect(refreshBtn);

      // Outer/inner border
      ctx.strokeStyle = '#444';
      ctx.lineWidth = 4 * scale;
      ctx.stroke();
      ctx.strokeStyle = '#666';
      ctx.lineWidth = 2 * scale;
      ctx.stroke();

      // Draw edges
      ctx.strokeStyle = GRAY;
      ctx.lineWidth = 2 * scale;
      edges.forEach(([i1, i2]) => {
        const { x: x1, y: y1 } = nodes[i1];
        const { x: x2, y: y2 } = nodes[i2];
        ctx.beginPath();
        ctx.moveTo(x1 * scale + offsetX, y1 * scale + offsetY);
        ctx.lineTo(x2 * scale + offsetX, y2 * scale + offsetY);
        ctx.stroke();
      });

      // Draw arrows
      arrows.forEach(({ from, to, color }) => {
        const { x: x1, y: y1 } = nodes[from];
        const { x: x2, y: y2 } = nodes[to];
        const dx = x2 - x1;
        const dy = y2 - y1;
        const length = Math.sqrt(dx * dx + dy * dy);
        const offset = 5 * scale;
        const offsetXperp = -dy / length * offset;
        const offsetYperp = dx / length * offset;
        const direction = color === 'blue' ? 1 : -1;
        const sx1 = x1 * scale + direction * offsetXperp + offsetX;
        const sy1 = y1 * scale + direction * offsetYperp + offsetY;
        const sx2 = x2 * scale + direction * offsetXperp + offsetX;
        const sy2 = y2 * scale + direction * offsetYperp + offsetY;
        ctx.strokeStyle = color;
        ctx.lineWidth = 2 * scale;
        ctx.beginPath();
        ctx.moveTo(sx1, sy1);
        ctx.lineTo(sx2, sy2);
        ctx.stroke();
        // Arrowhead
        const midX = (sx1 + sx2) / 2;
        const midY = (sy1 + sy2) / 2;
        const angle = Math.atan2(sy2 - sy1, sx2 - sx1);
        const arrowLength = 10 * scale;
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

      // Draw nodes
      nodes.forEach(({ x, y, clicks }, i) => {
        const drawX = x * scale + offsetX;
        const drawY = y * scale + offsetY;
        if (clicks === 0) ctx.fillStyle = LIGHT_GRAY;
        else if (clicks === 1) ctx.fillStyle = '#ff9999';
        else ctx.fillStyle = '#cc6666';
        ctx.beginPath();
        ctx.arc(drawX, drawY, RADIUS * scale, 0, Math.PI * 2);
        ctx.fill();
        ctx.strokeStyle = BLACK;
        ctx.lineWidth = 2 * scale;
        ctx.beginPath();
        ctx.arc(drawX, drawY, RADIUS * scale, 0, Math.PI * 2);
        ctx.stroke();
      });
    }, [nodes, arrows, refreshKey]);

    // --- FULL INTERACTIVE CLICK HANDLER WITH ARROWS (main copy.js logic) ---
    useEffect(() => {
      const canvas = canvasRef.current;
      if (!canvas) return;
      const handleMouseDown = (e) => {
        const rect = canvas.getBoundingClientRect();
        const scale = Math.min(PANEL_WIDTH / BASE_WIDTH, PANEL_HEIGHT / BASE_HEIGHT);
        const offsetX = (PANEL_WIDTH - BASE_WIDTH * scale) / 2;
        const offsetY = (PANEL_HEIGHT - BASE_HEIGHT * scale) / 2;
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
          const drawX = node.x * scale + offsetX;
          const drawY = node.y * scale + offsetY;
          const dx = drawX - x;
          const dy = drawY - y;
          if (Math.sqrt(dx * dx + dy * dy) <= RADIUS * scale) {
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
      return () => canvas.removeEventListener('mousedown', handleMouseDown);
    }, [refreshBtnRect, initialNodes, edges, maxClicks, RADIUS]);

    return (
      <canvas
        ref={canvasRef}
        width={PANEL_WIDTH}
        height={PANEL_HEIGHT}
        style={{ display: 'block', margin: 'auto', backgroundColor: '#1c1c1c', border: '1px solid #888', width: PANEL_WIDTH + 'px', height: PANEL_HEIGHT + 'px' }}
      />
    );
  };

  return (
    <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', width: '100%' }}>
      <div style={{ marginBottom: '1rem', width: panel_width }}>
        <input
          type="text"
          value={value}
          onChange={e => onChange(e.target.value)}
          placeholder="Filter sequences..."
          style={{ width: '100%', padding: '8px', borderRadius: '8px', border: '1px solid #444', background: '#222', color: '#fff' }}
        />
      </div>
      <Panel>
        <GameScaffoldPanel />
      </Panel>
    </div>
  );
};

export default SequenceFilter;
