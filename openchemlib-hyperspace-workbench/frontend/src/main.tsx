import {createRoot} from 'react-dom/client';
import App from './App';
import {AppearanceProvider} from './Appearance';
import './style.css';
import './precision.css';
createRoot(document.getElementById('root')!).render(<AppearanceProvider><App/></AppearanceProvider>);
