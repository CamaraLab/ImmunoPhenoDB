import React from 'react'
import ReactDOM from 'react-dom/client'
import App from './App.tsx'
import './index.css'
import 'bootstrap/dist/css/bootstrap.min.css';
import ReactGA from 'react-ga4'

// Initialize Google Analytics with your Measurement ID
ReactGA.initialize('G-XE1L8MNK9Q');

// Send a pageview when the app initially loads
ReactGA.send({ hitType: "pageview", page: window.location.pathname, title: "Initial Load" });

ReactDOM.createRoot(document.getElementById('root')!).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>,
)
