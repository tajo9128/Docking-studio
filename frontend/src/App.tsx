import { BrowserRouter, Routes, Route } from 'react-router-dom'
import { Layout } from '@/components/Layout'
import { Dashboard } from '@/pages/Dashboard'
import { Docking } from '@/pages/Docking'
import { JobQueue } from '@/pages/JobQueue'
import { Results } from '@/pages/Results'
import { AIAssistant } from '@/pages/AIAssistant'
import { Settings } from '@/pages/Settings'
import { MoleculeDynamics } from '@/pages/MoleculeDynamics'

export default function App() {
  return (
    <BrowserRouter>
      <Routes>
        <Route element={<Layout />}>
          <Route path="/" element={<Dashboard />} />
          <Route path="/docking" element={<Docking />} />
          <Route path="/md" element={<MoleculeDynamics />} />
          <Route path="/jobs" element={<JobQueue />} />
          <Route path="/results" element={<Results />} />
          <Route path="/ai" element={<AIAssistant />} />
          <Route path="/settings" element={<Settings />} />
        </Route>
      </Routes>
    </BrowserRouter>
  )
}
