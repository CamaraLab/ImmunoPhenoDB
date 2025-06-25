import { useState } from 'react'
import reactLogo from './assets/react.svg'
import viteLogo from '/vite.svg'
import './App.css'
import PlotlyLoad from './components/PlotlyLoad'
import MainPage from './components/MainPage'
import NavBar from './components/NavBar'

function App() {

  return (
    <>
      <NavBar />
      <div>
        <MainPage />
      </div>
    </>
  )
}

export default App
