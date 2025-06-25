import React from 'react'

import './css/plotly.css'
const PlotlyLoad = () => {
  return (
    <div>
        <div className="plotlybars-wrapper">
        <div className="plotlybars">
            <div className="plotlybars-bar b1"></div>
            <div className="plotlybars-bar b2"></div>
            <div className="plotlybars-bar b3"></div>
            <div className="plotlybars-bar b4"></div>
            <div className="plotlybars-bar b5"></div>
            <div className="plotlybars-bar b6"></div>
            <div className="plotlybars-bar b7"></div>
        </div>
        <div className="plotlybars-text">
            Loading plot...
        </div>
        </div>
    </div>
  )
}

export default PlotlyLoad