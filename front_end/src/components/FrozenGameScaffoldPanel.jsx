import React, { useRef, useEffect } from 'react';

const DEFAULT_NODES = [
  { x: 300, y: 300 },
  { x: 250, y: 250 }, { x: 350, y: 250 },
  { x: 200, y: 200 }, { x: 300, y: 200 }, { x: 400, y: 200 },
  { x: 150, y: 150 }, { x: 250, y: 150 }, { x: 350, y: 150 }, { x: 450, y: 150 },
  { x: 100, y: 100 }, { x: 200, y: 100 }, { x: 300, y: 100 }, { x: 400, y: 100 }, { x: 500, y: 100 }
];
const DEFAULT_EDGES = [
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

const PANEL_WIDTH = 570;
const PANEL_HEIGHT = 370;
const RADIUS = 20;
const LIGHT_GRAY = '#dddddd';
const GRAY = '#aaaaaa';
const BLACK = '#000000';

const FrozenGameScaffoldPanel = ({ initialNodes = DEFAULT_NODES, edges = DEFAULT_EDGES, panelWidth, panelHeight }) => {
  const canvasRef = useRef(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const width = panelWidth || PANEL_WIDTH;
    const height = panelHeight || PANEL_HEIGHT;
    canvas.width = width;
    canvas.height = height;
    const ctx = canvas.getContext('2d');
    ctx.clearRect(0, 0, width, height);
    // Calculate scaling factor to fit scaffold into the panel
    const scaleX = width / PANEL_WIDTH;
    const scaleY = height / PANEL_HEIGHT;
    const scale = Math.min(scaleX, scaleY);
    ctx.save();
    ctx.scale(scale, scale);
    // Center the scaffold
    const offsetX = (width / scale - PANEL_WIDTH) / 2;
    const offsetY = (height / scale - PANEL_HEIGHT) / 2;
    ctx.translate(offsetX, offsetY);
    // Draw background box
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
    // Draw edges
    ctx.strokeStyle = GRAY;
    ctx.lineWidth = 2;
    edges.forEach(([i1, i2]) => {
      const { x: x1, y: y1 } = initialNodes[i1];
      const { x: x2, y: y2 } = initialNodes[i2];
      ctx.beginPath();
      ctx.moveTo(x1, y1);
      ctx.lineTo(x2, y2);
      ctx.stroke();
    });
    // Draw nodes
    initialNodes.forEach(({ x, y }, i) => {
      ctx.fillStyle = LIGHT_GRAY;
      ctx.beginPath();
      ctx.arc(x, y, RADIUS, 0, Math.PI * 2);
      ctx.fill();
      ctx.strokeStyle = BLACK;
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(x, y, RADIUS, 0, Math.PI * 2);
      ctx.stroke();
    });
    ctx.restore();
  }, [initialNodes, edges, panelWidth, panelHeight]);

  return (
    <div style={{ position: 'relative', width: '100%', height: '100%' }}>
      <canvas ref={canvasRef} style={{ width: '100%', height: '100%', display: 'block' }} />
    </div>
  );
};

export default FrozenGameScaffoldPanel;
