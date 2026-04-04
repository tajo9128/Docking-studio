import { useState, useEffect } from 'react'
import { useTheme } from '@/contexts/ThemeContext'
import { getLocale, setLocale, locales, type Locale } from '@/i18n'

export function Header() {
  const { theme, toggleTheme } = useTheme()
  const isDark = theme === 'dark'
  const [language, setLanguageState] = useState<Locale>(getLocale)
  const [collapsed, setCollapsed] = useState(true)

  useEffect(() => {
    const handleStorage = () => setLanguageState(getLocale())
    window.addEventListener('storage', handleStorage)
    return () => window.removeEventListener('storage', handleStorage)
  }, [])

  const handleLanguageChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    const newLang = e.target.value as Locale
    setLocale(newLang)
    setLanguageState(newLang)
  }
  
  return (
    <header className={`shadow-sm transition-all ${
      isDark ? 'bg-gray-800 text-white' : 'bg-white text-gray-800 border-b border-gray-200'
    } ${collapsed ? 'h-10' : 'h-14'}`}>
      <div className="flex items-center justify-between px-3 h-full">
        <div className="flex items-center gap-3">
          <button
            onClick={() => setCollapsed(!collapsed)}
            className={`p-1 rounded transition-colors ${
              isDark ? 'hover:bg-gray-700 text-gray-400' : 'hover:bg-gray-100 text-gray-500'
            }`}
            title={collapsed ? 'Expand header' : 'Collapse header'}
          >
            <svg className={`w-4 h-4 transition-transform ${collapsed ? '' : 'rotate-180'}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
            </svg>
          </button>
          <a href="/" className="flex items-center gap-1">
            <span className="text-lg font-bold text-blue-500">BioDockify</span>
            <span className={`text-lg font-bold ${isDark ? 'text-white' : 'text-gray-900'}`}>Studio</span>
          </a>
          {!collapsed && (
            <span className={`hidden lg:block text-xs ${isDark ? 'text-gray-500' : 'text-gray-400'}`}>
              AI Powered Autonomous Drug Discovery Platform
            </span>
          )}
        </div>
        
        <div className={`flex items-center gap-2 transition-opacity ${collapsed ? 'opacity-0 pointer-events-none absolute' : 'opacity-100'}`}>
          <select
            value={language}
            onChange={handleLanguageChange}
            className={`px-2 py-1 text-sm border rounded transition-colors ${
              isDark
                ? 'bg-gray-700 border-gray-600 text-white'
                : 'bg-white border-gray-300 text-gray-900'
            }`}
            aria-label="Select language"
          >
            {locales.map((lang) => (
              <option key={lang.code} value={lang.code}>{lang.name}</option>
            ))}
          </select>

          <button
            onClick={toggleTheme}
            className={`p-1.5 rounded-lg transition-colors ${
              isDark
                ? 'hover:bg-gray-700 text-yellow-400'
                : 'hover:bg-gray-100 text-gray-600'
            }`}
            title={isDark ? 'Switch to light mode' : 'Switch to dark mode'}
          >
            {isDark ? (
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 3v1m0 16v1m9-9h-1M4 12H3m15.364 6.364l-.707-.707M6.343 6.343l-.707-.707m12.728 0l-.707.707M6.343 17.657l-.707.707M16 12a4 4 0 11-8 0 4 4 0 018 0z" />
              </svg>
            ) : (
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M20.354 15.354A9 9 0 018.646 3.646 9.003 9.003 0 0012 21a9.003 9.003 0 008.354-5.646z" />
              </svg>
            )}
          </button>
          
          <a
            href="/settings"
            className={`p-1.5 rounded-lg transition-colors ${
              isDark
                ? 'hover:bg-gray-700 text-gray-300'
                : 'hover:bg-gray-100 text-gray-600'
            }`}
            title="Settings"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10.325 4.317c.426-1.756 2.924-1.756 3.35 0a1.724 1.724 0 002.573 1.066c1.543-.94 3.31.826 2.37 2.37a1.724 1.724 0 001.065 2.572c1.756.426 1.756 2.924 0 3.35a1.724 1.724 0 00-1.066 2.573c.94 1.543-.826 3.31-2.37 2.37a1.724 1.724 0 00-2.572 1.065c-.426 1.756-2.924 1.756-3.35 0a1.724 1.724 0 00-2.573-1.066c-1.543.94-3.31-.826-2.37-2.37a1.724 1.724 0 00-1.065-2.572c-1.756-.426-1.756-2.924 0-3.35a1.724 1.724 0 001.066-2.573c-.94-1.543.826-3.31 2.37-2.37.996.608 2.296.07 2.572-1.065z" />
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
            </svg>
          </a>
        </div>
      </div>
    </header>
  )
}
