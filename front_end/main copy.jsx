import React, { useEffect, useRef } from 'react';
import './styles.css';

const GraphCanvas = () => {
  const canvasRef = useRef(null);
  const arrows = useRef([]);
  const RADIUS = 20;
  const LIGHT_GRAY = "#dddddd";
  const GRAY = "#aaaaaa";
  const BLACK = "#000000";

  const nodes = useRef([
    { x: 400, y: 300, clicks: 0 }, // 0 (root)
    { x: 350, y: 250, clicks: 0 }, { x: 450, y: 250, clicks: 0 },
    { x: 300, y: 200, clicks: 0 }, { x: 400, y: 200, clicks: 0 }, { x: 500, y: 200, clicks: 0 },
    { x: 250, y: 150, clicks: 0 }, { x: 350, y: 150, clicks: 0 }, { x: 450, y: 150, clicks: 0 }, { x: 550, y: 150, clicks: 0 },
    { x: 200, y: 100, clicks: 0 }, { x: 300, y: 100, clicks: 0 }, { x: 400, y: 100, clicks: 0 }, { x: 500, y: 100, clicks: 0 }, { x: 600, y: 100, clicks: 0 }
  ]);

  const edges = [
    [0, 1], [0, 2], [1, 3], [1, 4], [2, 4], [2, 5],
    [3, 6], [3, 7], [4, 7], [4, 8], [5, 8], [5, 9],
    [6, 10], [6, 11], [7, 11], [7, 12], [8, 12], [8, 13], [9, 13], [9, 14]
  ];

  const clickedNodes = useRef(new Set([0]));
  const lastNodeIndex = useRef(0);
  const maxClicks = 5;
  nodes.current[0].clicks = 1;

  useEffect(() => {
    const canvas = canvasRef.current;
    const ctx = canvas.getContext('2d');

    const draw = () => {
      const BASE_WIDTH = 900, BASE_HEIGHT = 600;
      canvas.width = BASE_WIDTH;
      canvas.height = BASE_HEIGHT;
      ctx.clearRect(0, 0, canvas.width, canvas.height);

      const xs = nodes.current.map(n => n.x);
      const ys = nodes.current.map(n => n.y);
      const minX = Math.min(...xs) - RADIUS - 30;
      const minY = Math.min(...ys) - RADIUS - 30;
      const offsetX = -minX + 100;
      const offsetY = -minY + 20;

      ctx.fillStyle = "#101010";
      ctx.fillRect(0, 0, canvas.width, canvas.height);

      ctx.strokeStyle = GRAY;
      ctx.lineWidth = 2;
      edges.forEach(([i1, i2]) => {
        const { x: x1, y: y1 } = nodes.current[i1];
        const { x: x2, y: y2 } = nodes.current[i2];
        ctx.beginPath();
        ctx.moveTo(x1 + offsetX, y1 + offsetY);
        ctx.lineTo(x2 + offsetX, y2 + offsetY);
        ctx.stroke();
      });

      arrows.current.forEach(({ from, to, color }) => {
        const { x: x1, y: y1 } = nodes.current[from];
        const { x: x2, y: y2 } = nodes.current[to];
        const dx = x2 - x1, dy = y2 - y1;
        const length = Math.sqrt(dx * dx + dy * dy);
        const offset = 5;
        const offsetXperp = -dy / length * offset;
        const offsetYperp = dx / length * offset;
        const direction = (color === "blue") ? 1 : -1;
        const sx1 = x1 + direction * offsetXperp + offsetX;
        const sy1 = y1 + direction * offsetYperp + offsetY;
        const sx2 = x2 + direction * offsetXperp + offsetX;
        const sy2 = y2 + direction * offsetYperp + offsetY;

        ctx.strokeStyle = color;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(sx1, sy1);
        ctx.lineTo(sx2, sy2);
        ctx.stroke();
      });

      nodes.current.forEach(({ x, y, clicks }) => {
        const drawX = x + offsetX, drawY = y + offsetY;
        ctx.fillStyle = clicks === 0 ? LIGHT_GRAY : (clicks === 1 ? "#ff9999" : "#cc6666");
        ctx.beginPath();
        ctx.arc(drawX, drawY, RADIUS, 0, Math.PI * 2);
        ctx.fill();
        ctx.strokeStyle = BLACK;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(drawX, drawY, RADIUS, 0, Math.PI * 2);
        ctx.stroke();
      });

      requestAnimationFrame(draw);
    };

    const handleClick = (e) => {
      const rect = canvas.getBoundingClientRect();
      const x = e.clientX - rect.left;
      const y = e.clientY - rect.top;
      const offsetX = 100, offsetY = 20;

      for (let i = 0; i < nodes.current.length; i++) {
        const node = nodes.current[i];
        const dx = node.x + offsetX - x;
        const dy = node.y + offsetY - y;
        if (Math.sqrt(dx * dx + dy * dy) <= RADIUS) {
          let isAllowed = false, proposedPrev = null;
          clickedNodes.current.forEach(prevIndex => {
            const isAdjacent = edges.some(([a, b]) => (a === prevIndex && b === i) || (b === prevIndex && a === i));
            const isAbove = node.y < nodes.current[prevIndex].y;
            const exists = arrows.current.some(({ from, to }) => from === prevIndex && to === i);
            if (isAdjacent && isAbove && !exists) {
              isAllowed = true;
              proposedPrev = prevIndex;
            }
          });

          if (!isAllowed) return;
          const totalClicks = nodes.current.reduce((sum, node) => sum + node.clicks, 0);
          if (totalClicks >= maxClicks) return;

          node.clicks += 1;
          clickedNodes.current.add(i);
          const from = proposedPrev ?? i;
          const to = i;
          if (from !== to) {
            const prevNode = nodes.current[from];
            const color = node.x < prevNode.x ? "red" : "blue";
            arrows.current.push({ from, to, color });
          }

          lastNodeIndex.current = i;
          break;
        }
      }
    };

    canvas.addEventListener("mousedown", handleClick);
    draw();
    return () => canvas.removeEventListener("mousedown", handleClick);
  }, []);

  return <canvas ref={canvasRef} style={{ width: '100%', height: '100%' }} />;
};

const Dashboard = () => {
  return (
    <div style={{
      display: 'grid',
      gridTemplateColumns: '2fr 1fr',
      gridTemplateRows: '2fr 1fr',
      gap: '10px',
      height: '100vh',
      backgroundColor: '#121212',
      padding: '10px',
      boxSizing: 'border-box',
      fontFamily: 'Inter, sans-serif'
    }}>
      <div style={{ backgroundColor: '#181818', borderRadius: '16px', overflow: 'hidden' }}>
        <GraphCanvas />
      </div>
      <div style={{ backgroundColor: '#181818', borderRadius: '16px' }}></div>
      <div style={{ backgroundColor: '#181818', borderRadius: '16px' }}></div>
      <div style={{ backgroundColor: '#181818', borderRadius: '16px' }}></div>
    </div>
  );
};

export default Dashboard;
