# 🚀 Quick Start Guide for Students

## Step 1: Start Docker Desktop

1. Open **Docker Desktop** on your computer
2. Wait for it to say "Docker is running"

## Step 2: Run the Backend

Open terminal/command prompt and run:

```bash
docker run -d -p 8000:8000 tajo9128/docking-studio:latest
```

## Step 3: Open in Browser

After starting the container, open your browser and go to:

### 📚 API Documentation (Recommended for learning)
➡️ **http://localhost:8000/docs**

This is the **Swagger UI** where you can:
- See all available API endpoints
- Test docking operations
- View request/response formats

### 🔧 Alternative: ReDoc
➡️ **http://localhost:8000/redoc**

### ✅ Health Check
➡️ **http://localhost:8000/health**

---

## Step 4: Run the Desktop App (Optional)

For the full PyQt6 desktop interface:

```bash
# Clone the repository
git clone https://github.com/tajo9128/Docking-studio.git
cd Docking-studio

# Install dependencies
pip install -r requirements.txt

# Run the app
python -m src.biodockify_main
```

The desktop app will connect to `http://localhost:8000` automatically.

---

## 🎯 Quick Reference

| Service | URL |
|---------|-----|
| Swagger API Docs | http://localhost:8000/docs |
| ReDoc | http://localhost:8000/redoc |
| Health Status | http://localhost:8000/health |
| Security Status | http://localhost:8000/security/status |

---

## ❓ Troubleshooting

**Port already in use?**
```bash
docker stop $(docker ps -q)
docker run -d -p 8000:8000 tajo9128/docking-studio:latest
```

**Check if running:**
```bash
docker ps
```

**View logs:**
```bash
docker logs $(docker ps -q --last-format)
```

---

**Happy Docking! 🧬**
