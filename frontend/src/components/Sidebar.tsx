import { NavLink } from 'react-router-dom'
import { clsx } from 'clsx'

const navItems = [
  {
    section: 'MAIN',
    items: [
      { path: '/', label: 'Dashboard', icon: '🏠' },
      { path: '/docking', label: 'New Docking', icon: '🧪' },
      { path: '/jobs', label: 'Job Queue', icon: '📋' },
    ],
  },
  {
    section: 'ANALYSIS',
    items: [
      { path: '/results', label: 'Results', icon: '📊' },
      { path: '/rmsd', label: 'RMSD Analysis', icon: '🎯' },
      { path: '/interactions', label: 'Interactions', icon: '🔗' },
      { path: '/pharmacophore', label: 'Pharmacophore', icon: '🧬' },
      { path: '/qsar', label: 'QSAR Modeling', icon: '🧮' },
    ],
  },
  {
    section: 'TOOLS',
    items: [
      { path: '/viewer', label: 'Molecule Viewer', icon: '👁' },
      { path: '/ai', label: 'AI Assistant', icon: '🤖' },
      { path: '/security', label: 'Security', icon: '🔒' },
    ],
  },
]

export function Sidebar() {
  return (
    <aside className="w-60 h-full bg-background-sidebar flex flex-col">
      {/* Logo */}
      <div className="px-5 py-6 bg-background-dark">
        <h1 className="text-white text-xl font-bold">Docking Studio</h1>
        <p className="text-gray-500 text-xs mt-0.5">Version 2.0.1</p>
      </div>

      {/* Navigation */}
      <nav className="flex-1 overflow-y-auto py-4">
        {navItems.map((section) => (
          <div key={section.section} className="mb-4">
            <p className="px-5 text-gray-500 text-xs font-semibold tracking-wider mb-2">
              {section.section}
            </p>
            {section.items.map((item) => (
              <NavLink
                key={item.path}
                to={item.path}
                end={item.path === '/'}
                className={({ isActive }) =>
                  clsx(
                    'flex items-center gap-3 px-5 py-2.5 text-sm font-medium transition-colors',
                    isActive
                      ? 'bg-primary text-white'
                      : 'text-gray-300 hover:bg-gray-800 hover:text-white'
                  )
                }
              >
                <span>{item.icon}</span>
                {item.label}
              </NavLink>
            ))}
          </div>
        ))}
      </nav>

      {/* Bottom section */}
      <div className="p-4 border-t border-gray-800 space-y-1">
        <NavLink
          to="/settings"
          className="flex items-center gap-3 px-3 py-2 text-gray-400 hover:text-white hover:bg-gray-800 rounded-lg transition-colors text-sm"
        >
          <span>⚙</span>
          AI Settings
        </NavLink>
        <NavLink
          to="/settings"
          className="flex items-center gap-3 px-3 py-2 text-gray-400 hover:text-white hover:bg-gray-800 rounded-lg transition-colors text-sm"
        >
          <span>🔧</span>
          Preferences
        </NavLink>
      </div>
    </aside>
  )
}
