let arrows = [];
const canvas = document.getElementById('graphCanvas');
const ctx = canvas.getContext('2d');

const RADIUS = 20;
const LIGHT_GRAY = "#dddddd";
const GRAY = "#aaaaaa";
const BLACK = "#000000";

let nodes = [
    { x: 400, y: 300, clicks: 0 }, // 0 (root node pre-colored pink)
    { x: 350, y: 250, clicks: 0 }, { x: 450, y: 250, clicks: 0 },
    { x: 300, y: 200, clicks: 0 }, { x: 400, y: 200, clicks: 0 }, { x: 500, y: 200, clicks: 0 },
    { x: 250, y: 150, clicks: 0 }, { x: 350, y: 150, clicks: 0 }, { x: 450, y: 150, clicks: 0 }, { x: 550, y: 150, clicks: 0 },
    { x: 200, y: 100, clicks: 0 }, { x: 300, y: 100, clicks: 0 }, { x: 400, y: 100, clicks: 0 }, { x: 500, y: 100, clicks: 0 }, { x: 600, y: 100, clicks: 0 }
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
    ctx.clearRect(0, 0, canvas.width, canvas.height);

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

    // Draw directional arrows
    // These arrows represent a directed linkage from one ubiquitin to another,
    // e.g. simulating a chain where one ubiquitin is bound to another at a specific site.
    arrows.forEach(({ from, to, color }) => {
        const { x: x1, y: y1 } = nodes[from];
        const { x: x2, y: y2 } = nodes[to];

        // Calculate perpendicular offset
        const dx = x2 - x1;
        const dy = y2 - y1;
        const length = Math.sqrt(dx * dx + dy * dy);
        const offset = 5;
        const offsetX = -dy / length * offset;
        const offsetY = dx / length * offset;
        const direction = (color === "blue") ? 1 : -1;

        const sx1 = x1 + direction * offsetX;
        const sy1 = y1 + direction * offsetY;
        const sx2 = x2 + direction * offsetX;
        const sy2 = y2 + direction * offsetY;

        // Draw the line
        ctx.strokeStyle = color;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(sx1, sy1);
        ctx.lineTo(sx2, sy2);
        ctx.stroke();

        // Draw arrowhead at midpoint
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
    nodes.forEach(({ x, y, clicks }) => {
        if (clicks === 0) ctx.fillStyle = LIGHT_GRAY;
        else if (clicks === 1) ctx.fillStyle = "#ff9999";
        else ctx.fillStyle = "#cc6666";
        ctx.beginPath();
        ctx.arc(x, y, RADIUS, 0, Math.PI * 2);
        ctx.fill();

        ctx.strokeStyle = BLACK;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(x, y, RADIUS, 0, Math.PI * 2);
        ctx.stroke();
    });

    requestAnimationFrame(draw);
}

// Handle clicks
canvas.addEventListener('mousedown', (e) => {
    const rect = canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;

    for (let i = 0; i < nodes.length; i++) {
        const node = nodes[i];
        const dx = node.x - x;
        const dy = node.y - y;

        if (Math.sqrt(dx * dx + dy * dy) <= RADIUS) {
            let isAllowed = false;
            let proposedPrev = null;
            clickedNodes.forEach(prevIndex => {
                const isAdjacent = edges.some(([a, b]) =>
                    (a === prevIndex && b === i) || (b === prevIndex && a === i)
                );
                const isAbove = nodes[i].y < nodes[prevIndex].y;
                const lineExists = arrows.some(({ from, to }) =>
                    (from === prevIndex && to === i)
                );
                if (isAdjacent && isAbove && !lineExists) {
                    isAllowed = true;
                    proposedPrev = prevIndex;
                }
            });

            // Allow clicking the root node (i === 0) if it meets adjacency and direction rules
            if (!isAllowed) continue;

            const totalClicks = nodes.reduce((sum, node) => sum + node.clicks, 0);
            if (totalClicks >= maxClicks) return;

            node.clicks = (node.clicks || 0) + 1;
            if (!clickedNodes.has(i)) {
                clickedNodes.add(i);
            }

            // Use proposedPrev for lastNodeIndex
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
});

draw();

// Export linkage data functionality
function exportLinkages() {
    const linkages = arrows.map(({ from, to, color }) => ({
        from,
        to,
        linkage: color === 'red' ? 'K48' : 'K63'
    }));

    const blob = new Blob([JSON.stringify(linkages, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = "linkages_output.json";
    a.click();
    URL.revokeObjectURL(url);
}

const exportButton = document.createElement('button');
exportButton.innerText = "Export Linkages";
exportButton.style.position = "absolute";
exportButton.style.top = "10px";
exportButton.style.left = "10px";
exportButton.addEventListener("click", exportLinkages);
document.body.appendChild(exportButton);

