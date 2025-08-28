function simulateClicksFromJson(json, nodes, edges) {
  let clicked = new Set();
  let arrows = [];
  let mapping = {};
  let smallNodes = [];

  function recurse(ub, parentIdx, used, originIdx = parentIdx) {
    nodes[parentIdx] = { ...nodes[parentIdx], clicks: 1, chain_number: ub.chain_number };
    mapping[ub.chain_number] = parentIdx;
    clicked.add(parentIdx);

    for (const site of ub.branching_sites || []) {
      if (site.children && typeof site.children === 'object') {
        const childIdx = findNextAllowedNode(parentIdx, used, site.site_name);
        if (childIdx !== null && childIdx < nodes.length) {
          let color = 'gray';
          if (site.site_name === 'K48') color = 'red';
          else if (site.site_name === 'K63') color = 'blue';

          arrows.push({ from: parentIdx, to: childIdx, color, highlighted: true });
          used.add(childIdx);
          recurse(site.children, childIdx, used, originIdx);
        }
      }
    }

    const originalIdx = mapping[ub.chain_number];
    if (originalIdx !== undefined && originalIdx !== parentIdx) {
      const originalNode = nodes[originalIdx];
      const availableNodes = nodes.filter((node, idx) =>
        !used.has(idx) &&
        node.y === originalNode.y - 50 &&
        Math.abs(node.x - originalNode.x) <= 100
      );
      for (const nextNode of availableNodes) {
        const nextNodeIdx = nodes.indexOf(nextNode);
        if (!used.has(nextNodeIdx)) {
          arrows.push({ from: originalIdx, to: nextNodeIdx, color: 'gray', highlighted: true });
          used.add(nextNodeIdx);
          break;
        }
      }
    }

    for (const site of ub.branching_sites || []) {
      if (site.children === 'SMAC' || site.children === 'ABOC') {
        const colorState = site.children === 'SMAC' ? 2 : 1;
        const angle = site.site_name === 'K48' ? -Math.PI / 4 : Math.PI / 4;
        const dx = Math.sin(angle) * 25;
        const dy = -Math.cos(angle) * 25;
        const node = nodes[parentIdx];
        smallNodes.push({
          x: node.x + dx,
          y: node.y + dy,
          key: `protect-${node.chain_number}-${site.site_name}`,
          state: colorState
        });
      }
    }
  }

  function findNextAllowedNode(parentIdx, used, siteName) {
    // K63 moves rightward and up, K48 moves leftward and up
    const direction = siteName === 'K63' ? 1 : -1;
    const targetY = nodes[parentIdx].y - 50;
    for (let i = 0; i < nodes.length; i++) {
      if (used.has(i)) continue;

      const dx = nodes[i].x - nodes[parentIdx].x;
      const dy = nodes[i].y - nodes[parentIdx].y;
      const isCorrectDirection = direction === -1 ? dx < 0 : dx > 0;

      const isDirectlyAbove = nodes[i].y === targetY && Math.abs(dx) <= 100;
      const isAdjacent = edges.some(([a, b]) => (a === parentIdx && b === i) || (b === parentIdx && a === i));

      if (isDirectlyAbove && isAdjacent && isCorrectDirection) {
        return i;
      }
    }
    return null;
  }

  let used = new Set([0]);
  recurse(json, 0, used);
  return { nodes, arrows, mapping, smallNodes };
}

export { simulateClicksFromJson };