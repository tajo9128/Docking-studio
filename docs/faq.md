# Frequently Asked Questions (FAQ)

**Q: Is BioDockify free?**
A: Yes, BioDockify is open-source software released under the Apache 2.0 License.

**Q: Can I run this without Docker?**
A: No. BioDockify relies on Docker containers to ensure consistent scientific environments across all operating systems.

**Q: Why does the first job take so long?**
A: The first run downloads the docking engine image (~500MB). Subsequent runs will be much faster.

**Q: What docking engine is used?**
A: We use AutoDock Vina 1.2.3 for the core simulation, with RDKit for preprocessing and ODDT for scoring.

**Q: How do I report a bug?**
A: Please open an issue on our [GitHub Repository](https://github.com/tajo9128/Docking-studio).
