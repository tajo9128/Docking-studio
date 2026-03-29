import { BrowserRouter, Routes, Route } from 'react-router-dom'
import { Layout } from '@/components/Layout'
import { Dashboard } from '@/pages/Dashboard'
import { Docking } from '@/pages/Docking'
import { JobQueue } from '@/pages/JobQueue'
import { Results } from '@/pages/Results'
import { RMSDAnalysis } from '@/pages/RMSDAnalysis'
import { Interactions } from '@/pages/Interactions'
import { Viewer } from '@/pages/Viewer'
import { AIAssistant } from '@/pages/AIAssistant'
import { Security } from '@/pages/Security'
import { Settings } from '@/pages/Settings'
import { Pharmacophore } from '@/pages/Pharmacophore'
import { QSARModeling } from '@/pages/QSARModeling'

export default function App() {
  return (
    <BrowserRouter>
      <Routes>
        <Route element={<Layout />}>
          <Route path="/" element={<Dashboard />} />
          <Route path="/docking" element={<Docking />} />
          <Route path="/jobs" element={<JobQueue />} />
          <Route path="/results" element={<Results />} />
          <Route path="/rmsd" element={<RMSDAnalysis />} />
          <Route path="/interactions" element={<Interactions />} />
          <Route path="/pharmacophore" element={<Pharmacophore />} />
          <Route path="/qsar" element={<QSARModeling />} />
          <Route path="/viewer" element={<Viewer />} />
          <Route path="/ai" element={<AIAssistant />} />
          <Route path="/security" element={<Security />} />
          <Route path="/settings" element={<Settings />} />
        </Route>
      </Routes>
    </BrowserRouter>
  )
}
