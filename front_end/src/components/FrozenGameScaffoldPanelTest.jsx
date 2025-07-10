import React, { useRef, useEffect, useState } from 'react';

const DEFAULT_NODES = [
  { x: 300, y: 300, clicks: 0 },
  { x: 250, y: 250, clicks: 0 }, { x: 350, y: 250, clicks: 0 },
  { x: 200, y: 200, clicks: 0 }, { x: 300, y: 200, clicks: 0 }, { x: 400, y: 200, clicks: 0 },
  { x: 150, y: 150, clicks: 0 }, { x: 250, y: 150, clicks: 0 }, { x: 350, y: 150, clicks: 0 }, { x: 450, y: 150, clicks: 0 },
  { x: 100, y: 100, clicks: 0 }, { x: 200, y: 100, clicks: 0 }, { x: 300, y: 100, clicks: 0 }, { x: 400, y: 100, clicks: 0 }, { x: 500, y: 100, clicks: 0 }
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

const FrozenGameScaffoldPanelTest = ({ initialNodes = DEFAULT_NODES, edges = DEFAULT_EDGES, panelWidth, panelHeight, arrows: arrowsProp, frozen = false }) => {
  const canvasRef = useRef(null);
  const [nodes, setNodes] = useState(() => {
    const n = initialNodes.map(node => ({ ...node }));
    n[0].clicks = 1;
    return n;
  });
  const [arrows, setArrows] = useState(arrowsProp || []);
  const [clickedNodes, setClickedNodes] = useState(new Set([0]));
  const [lastNodeIndex, setLastNodeIndex] = useState(0);
  const maxClicks = 5;

  // Small node state: array of {x, y, state} where state: 0=none, 1=red, 2=blue
  const [smallNodes, setSmallNodes] = useState([]);
  const smallNodeRadius = RADIUS / 4;

  // refs for nodes, arrows, clickedNodes, smallNodes
  const nodesRef = useRef(nodes);
  const arrowsRef = useRef(arrows);
  const clickedNodesRef = useRef(clickedNodes);
  const smallNodesRef = useRef(smallNodes);
  useEffect(() => { nodesRef.current = nodes; }, [nodes]);
  useEffect(() => { arrowsRef.current = arrows; }, [arrows]);
  useEffect(() => { clickedNodesRef.current = clickedNodes; }, [clickedNodes]);
  useEffect(() => { smallNodesRef.current = smallNodes; }, [smallNodes]);

  // Helper to generate all small node positions and keep their state
  useEffect(() => {
    // Collect all small node positions (with unique keys)
    let smalls = [];
    // 1. At the lower side of each upward edge
    edges.forEach(([i1, i2]) => {
      // Only add if the lower node has been clicked
      if (nodes[i1].clicks > 0) {
        const { x: x1, y: y1 } = nodes[i1];
        const { x: x2, y: y2 } = nodes[i2];
        if (y2 < y1) {
          const dx = x2 - x1;
          const dy = y2 - y1;
          const len = Math.sqrt(dx * dx + dy * dy);
          const t = (RADIUS + smallNodeRadius) / len;
          const cx = x1 + dx * t;
          const cy = y1 + dy * t;
          smalls.push({ x: cx, y: cy, key: `edge-${i1}-${i2}` });
        }
      }
    });
    // 2. On the outermost nodes (nodes with only one edge)
    nodes.forEach((node, idx) => {
      const connectedEdges = edges.filter(([a, b]) => a === idx || b === idx);
      if (connectedEdges.length === 1 && node.clicks > 0) {
        const edge = connectedEdges[0];
        const otherIdx = edge[0] === idx ? edge[1] : edge[0];
        let dx = node.x - nodes[otherIdx].x;
        let dy = node.y - nodes[otherIdx].y;
        const len = Math.sqrt(dx * dx + dy * dy);
        dx = dx / len;
        dy = dy / len;
        const cx = node.x + dx * (RADIUS + smallNodeRadius);
        const cy = node.y + dy * (RADIUS + smallNodeRadius);
        smalls.push({ x: cx, y: cy, key: `outer-${idx}` });
      }
    });
    // 3. For every node, two small nodes at ±45° from vertical
    nodes.forEach((node, idx) => {
      if (node.clicks > 0) {
        const angles = [-Math.PI / 4, Math.PI / 4];
        angles.forEach((angle, j) => {
          const dx = Math.sin(angle) * (RADIUS + smallNodeRadius);
          const dy = -Math.cos(angle) * (RADIUS + smallNodeRadius);
          const cx = node.x + dx;
          const cy = node.y + dy;
          smalls.push({ x: cx, y: cy, key: `angle-${idx}-${j}` });
        });
      }
    });
    // Preserve state for existing small nodes
    setSmallNodes(prev => {
      const stateMap = {};
      prev.forEach(sn => { stateMap[sn.key] = sn.state || 0; });
      return smalls.map(sn => ({ ...sn, state: stateMap[sn.key] || 0 }));
    });
  }, [nodes, edges]);

  // Draw the scaffold
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const width = panelWidth || PANEL_WIDTH;
    const height = panelHeight || PANEL_HEIGHT;
    canvas.width = width;
    canvas.height = height;
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    // Draw background box
    const boxX = 0;
    const boxY = 0;
    const boxW = width;
    const boxH = height;
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
    // Outer/inner border
    ctx.strokeStyle = '#444';
    ctx.lineWidth = 4;
    ctx.stroke();
    ctx.strokeStyle = '#666';
    ctx.lineWidth = 2;
    ctx.stroke();
    // Draw edges
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
    // Draw small nodes (with color by state)
    smallNodes.forEach(sn => {
      if (sn.state === 0) return; // skip if no state
      else if (sn.state === 1) ctx.fillStyle = 'rgba(255,100,100,0.7)'; // hazy red
      else if (sn.state === 2) ctx.fillStyle = 'rgba(100,100,255,0.7)'; // hazy blue
      else ctx.fillStyle = '#888';
      ctx.beginPath();
      ctx.arc(sn.x, sn.y, smallNodeRadius, 0, Math.PI * 2);
      ctx.fill();
      ctx.strokeStyle = BLACK;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.arc(sn.x, sn.y, smallNodeRadius, 0, Math.PI * 2);
      ctx.stroke();
    });
    // Draw arrows
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
    // Draw nodes
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
  }, [nodes, arrows, panelWidth, panelHeight, smallNodes]);

  // Click handler (disable if frozen)
  useEffect(() => {
    if (frozen) return; // Do not attach click handler if frozen
    const canvas = canvasRef.current;
    if (!canvas) return;
    const handleMouseDown = (e) => {
      const rect = canvas.getBoundingClientRect();
      const x = (e.clientX - rect.left);
      const y = (e.clientY - rect.top);
      // Check small nodes first
      for (let i = 0; i < smallNodesRef.current.length; i++) {
        const sn = smallNodesRef.current[i];
        const dx = sn.x - x;
        const dy = sn.y - y;
        if (Math.sqrt(dx * dx + dy * dy) <= smallNodeRadius) {
          setSmallNodes(prev => prev.map((n, idx) => idx === i ? { ...n, state: n.state === 0 ? 1 : n.state === 1 ? 2 : 0 } : n));
          return;
        }
      }
      // ...existing code for main node clicks...
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
  }, [nodes, arrows, clickedNodes, lastNodeIndex, smallNodes, frozen]);

  return <canvas ref={canvasRef} style={{ width: '100%', height: '100%', cursor: frozen ? 'not-allowed' : 'pointer' }} />;
};

export default FrozenGameScaffoldPanelTest;
