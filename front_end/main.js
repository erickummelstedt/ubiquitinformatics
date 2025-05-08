const panel_width = 600;
const panel_height = 400;
const canvas_width = panel_width-30;
const canvas_height = panel_height-30;

const dashboardStyle = {
  display: 'grid',
  gridTemplateColumns: 'auto auto',
  gridTemplateRows: 'auto auto',
  gap: '10px',
  padding: '10px',
  boxSizing: 'border-box'
};

const basePanelStyle = {
  backgroundColor: '#181818',
  padding: '10px',
  border: '1px solid',
  borderColor: '#444444',
  width: panel_width.toString() + 'px',
  height: panel_height.toString() + 'px',
  overflow: 'auto',
  boxSizing: 'border-box',
  color: '#E0E0E0',
  cursor: 'pointer',
  borderRadius: '16px'
};

const Panel = ({ children }) => {
  const [hover, setHover] = React.useState(false);
  const style = {
    ...basePanelStyle,
    borderColor: hover ? '#888888' : '#444444',
  };

  return React.createElement(
    'div',
    {
      style: style,
      onMouseEnter: () => setHover(true),
      onMouseLeave: () => setHover(false),
    },
    children
  );
};

const GraphPanel = () => {
  const canvasRef = React.useRef(null);

  React.useEffect(() => {
    let arrows = [];
    const canvas = canvasRef.current;
    const ctx = canvas.getContext('2d');

    const RADIUS = 20;
    const LIGHT_GRAY = "#dddddd";
    const GRAY = "#aaaaaa";
    const BLACK = "#000000";

    let nodes = [
        { x: 285, y: 300, clicks: 0 }, // 0 (root node pre-colored pink)
        { x: 235, y: 250, clicks: 0 }, { x: 335, y: 250, clicks: 0 },
        { x: 185, y: 200, clicks: 0 }, { x: 285, y: 200, clicks: 0 }, { x: 385, y: 200, clicks: 0 },
        { x: 135, y: 150, clicks: 0 }, { x: 235, y: 150, clicks: 0 }, { x: 335, y: 150, clicks: 0 }, { x: 435, y: 150, clicks: 0 },
        { x: 85, y: 100, clicks: 0 }, { x: 185, y: 100, clicks: 0 }, { x: 285, y: 100, clicks: 0 }, { x: 385, y: 100, clicks: 0 }, { x: 485, y: 100, clicks: 0 }
    ];

    let edges = [
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

    let lastNodeIndex = null;
    let maxClicks = 5;
    let clickedNodes = new Set();

    nodes[0].clicks = 1;
    clickedNodes.add(0);
    lastNodeIndex = 0;

    function draw() {
      const BASE_WIDTH = canvas_width;
      const BASE_HEIGHT = canvas_height;
      canvas.width = BASE_WIDTH;
      canvas.height = BASE_HEIGHT;
      ctx.clearRect(0, 0, canvas.width, canvas.height);

      const padding = 30;
      const xs = nodes.map(n => n.x);
      const ys = nodes.map(n => n.y);
      const minX = Math.min(...xs) - RADIUS - padding;
      const maxX = Math.max(...xs) + RADIUS + padding;
      const minY = Math.min(...ys) - RADIUS - padding;
      const maxY = Math.max(...ys) + RADIUS + padding;
      const offsetX = -minX + 100;
      const offsetY = -minY + 20;


      ctx.strokeStyle = GRAY;
      ctx.lineWidth = 2;
      edges.forEach(([i1, i2]) => {
        const { x: x1, y: y1 } = nodes[i1];
        const { x: x2, y: y2 } = nodes[i2];
        ctx.beginPath();
        ctx.moveTo(x1 + offsetX, y1 + offsetY);
        ctx.lineTo(x2 + offsetX, y2 + offsetY);
        ctx.stroke();
      });

      arrows.forEach(({ from, to, color }) => {
        const { x: x1, y: y1 } = nodes[from];
        const { x: x2, y: y2 } = nodes[to];
        const dx = x2 - x1;
        const dy = y2 - y1;
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

        const midX = (sx1 + sx2) / 2;
        const midY = (sy1 + sy2) / 2;
        const angle = Math.atan2(sy2 - sy1, sx2 - sx1);
        const arrowLength = 10;

        ctx.fillStyle = color;
        ctx.beginPath();
        ctx.moveTo(midX, midY);
        ctx.lineTo(midX - arrowLength * Math.cos(angle - Math.PI / 6), midY - arrowLength * Math.sin(angle - Math.PI / 6));
        ctx.lineTo(midX - arrowLength * Math.cos(angle + Math.PI / 6), midY - arrowLength * Math.sin(angle + Math.PI / 6));
        ctx.closePath();
        ctx.fill();
      });

      nodes.forEach(({ x, y, clicks }) => {
        const drawX = x + offsetX;
        const drawY = y + offsetY;
        if (clicks === 0) ctx.fillStyle = LIGHT_GRAY;
        else if (clicks === 1) ctx.fillStyle = "#ff9999";
        else ctx.fillStyle = "#cc6666";
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
    }

    const handleCanvasClick = (e) => {
      const rect = canvas.getBoundingClientRect();
      const x = e.clientX - rect.left;
      const y = e.clientY - rect.top;
      const xs = nodes.map(n => n.x);
      const ys = nodes.map(n => n.y);
      const minX = Math.min(...xs) - RADIUS - 30;
      const minY = Math.min(...ys) - RADIUS - 30;
      const offsetX = -minX + 100;
      const offsetY = -minY + 20;

      for (let i = 0; i < nodes.length; i++) {
        const node = nodes[i];
        const dx = node.x + offsetX - x;
        const dy = node.y + offsetY - y;

        if (Math.sqrt(dx * dx + dy * dy) <= RADIUS) {
          let isAllowed = false;
          let proposedPrev = null;
          clickedNodes.forEach(prevIndex => {
            const isAdjacent = edges.some(([a, b]) => (a === prevIndex && b === i) || (b === prevIndex && a === i));
            const isAbove = nodes[i].y < nodes[prevIndex].y;
            const lineExists = arrows.some(({ from, to }) => (from === prevIndex && to === i));
            if (isAdjacent && isAbove && !lineExists) {
              isAllowed = true;
              proposedPrev = prevIndex;
            }
          });

          if (!isAllowed) continue;

          const totalClicks = nodes.reduce((sum, node) => sum + node.clicks, 0);
          if (totalClicks >= maxClicks) return;

          node.clicks = (node.clicks || 0) + 1;
          if (!clickedNodes.has(i)) {
            clickedNodes.add(i);
          }

          lastNodeIndex = proposedPrev !== null ? proposedPrev : i;

          if (lastNodeIndex !== null && lastNodeIndex !== i) {
            const prevNode = nodes[lastNodeIndex];
            const newNode = nodes[i];
            const color = newNode.x < prevNode.x ? "red" : "blue";
            arrows.push({ from: lastNodeIndex, to: i, color });
          }

          lastNodeIndex = i;
          return;
        }
      }
    };

    canvas.addEventListener('mousedown', handleCanvasClick);

    draw();
  }, []);

  return React.createElement('canvas', {
    ref: canvasRef,
    width: 300,
    height: 200,
    style: { display: 'block', margin: 'auto', backgroundColor: '#1c1c1c', border: '1px solid #888', width : canvas_width.toString() + 'px', height: canvas_height.toString() + 'px',},
  });
};

const App = () =>
  React.createElement(
    'div',
    { style: dashboardStyle },
    React.createElement(
      Panel,
      null,
      React.createElement(GraphPanel)
    ),
    React.createElement(Panel, null, React.createElement('h2', null, 'Top Right Panel')),
    React.createElement(Panel, null, React.createElement('h2', null, 'Bottom Left Panel')),
    React.createElement(Panel, null, React.createElement('h2', null, 'Bottom Right Panel'))
  );

ReactDOM.createRoot(document.getElementById('root')).render(React.createElement(App));
