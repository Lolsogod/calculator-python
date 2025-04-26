# Как запускать

## Установка виртуального окружения и зависимостей

```bash
python3 -m venv venv
source venv/bin/activate  # Для Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Запуск приложения

```bash
uvicorn main:app --reload
```
