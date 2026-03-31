import { useState } from 'react'
import Dashboard from './pages/Dashboard'
import Docking from './pages/Docking'
import Results from './pages/Results'
import MD from './pages/MD'
import AI from './pages/AI'
import Settings from './pages/Settings'

type Page = 'dashboard' | 'docking' | 'results' | 'md' | 'ai' | 'settings'

export default function App() {
  const [currentPage, setCurrentPage] = useState<Page>('dashboard')

  const renderPage = () => {
    switch (currentPage) {
      case 'dashboard': return <Dashboard />
      case 'docking': return <Docking />
      case 'results': return <Results />
      case 'md': return <MD />
      case 'ai': return <AI />
      case 'settings': return <Settings />
      default: return <Dashboard />
    }
  }

  return (
    <div className="app-container">
      <header className="header">
        <h1>BioDockify</h1>
        <nav className="header-nav">
          <a
            href="#"
            className={currentPage === 'dashboard' ? 'active' : ''}
            onClick={e => { e.preventDefault(); setCurrentPage('dashboard') }}
          >
            Home
          </a>
          <a
            href="#"
            className={currentPage === 'docking' ? 'active' : ''}
            onClick={e => { e.preventDefault(); setCurrentPage('docking') }}
          >
            Docking
          </a>
          <a
            href="#"
            className={currentPage === 'results' ? 'active' : ''}
            onClick={e => { e.preventDefault(); setCurrentPage('results') }}
          >
            Results
          </a>
          <a
            href="#"
            className={currentPage === 'ai' ? 'active' : ''}
            onClick={e => { e.preventDefault(); setCurrentPage('ai') }}
          >
            AI Assistant
          </a>
        </nav>
      </header>

      <div className="main-content">
        <aside className="sidebar">
          <nav className="sidebar-nav">
            <a
              href="#"
              className={currentPage === 'dashboard' ? 'active' : ''}
              onClick={e => { e.preventDefault(); setCurrentPage('dashboard') }}
            >
              <span className="sidebar-icon">&#127968;</span>
              Dashboard
            </a>
            <a
              href="#"
              className={currentPage === 'docking' ? 'active' : ''}
              onClick={e => { e.preventDefault(); setCurrentPage('docking') }}
            >
              <span className="sidebar-icon">&#128300;</span>
              Docking
            </a>
            <a
              href="#"
              className={currentPage === 'md' ? 'active' : ''}
              onClick={e => { e.preventDefault(); setCurrentPage('md') }}
            >
              <span className="sidebar-icon">&#128200;</span>
              MD Simulation
            </a>
            <a
              href="#"
              className={currentPage === 'results' ? 'active' : ''}
              onClick={e => { e.preventDefault(); setCurrentPage('results') }}
            >
              <span className="sidebar-icon">&#128202;</span>
              Results
            </a>
            <a
              href="#"
              className={currentPage === 'ai' ? 'active' : ''}
              onClick={e => { e.preventDefault(); setCurrentPage('ai') }}
            >
              <span className="sidebar-icon">&#129302;</span>
              AI Assistant
            </a>
            <a
              href="#"
              className={currentPage === 'settings' ? 'active' : ''}
              onClick={e => { e.preventDefault(); setCurrentPage('settings') }}
            >
              <span className="sidebar-icon">&#9881;</span>
              Settings
            </a>
          </nav>
        </aside>

        <main className="page">
          {renderPage()}
        </main>
      </div>

      <footer className="status-bar">
        <div className="status-item">
          <span className="status-dot"></span>
          <span>Connected to BioDockify Backend</span>
        </div>
        <div className="status-item">
          <span>v2.3.0</span>
        </div>
      </footer>
    </div>
  )
}
