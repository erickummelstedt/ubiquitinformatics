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
  width: '600px',
  height: '400px',
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

const App = () =>
  React.createElement(
    'div',
    { style: dashboardStyle },
    React.createElement(
      Panel,
      null,
      React.createElement('h1', null, 'Hello from CDN React!'),
      React.createElement('p', null, 'Hello World')
    ),
    React.createElement(Panel, null, React.createElement('h2', null, 'Top Right Panel')),
    React.createElement(Panel, null, React.createElement('h2', null, 'Bottom Left Panel')),
    React.createElement(Panel, null, React.createElement('h2', null, 'Bottom Right Panel'))
  );

ReactDOM.createRoot(document.getElementById('root')).render(React.createElement(App));
