#!/usr/bin/env bash
set -euo pipefail

DB_NAME="${1:-labbench2pro}"
DB_USER="${2:-dev}"
DB_PORT="${3:-5433}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Creating database '$DB_NAME' on port $DB_PORT (if not exists)..."
createdb -h localhost -p "$DB_PORT" -U "$DB_USER" "$DB_NAME" 2>/dev/null || echo "Database already exists."

echo "Applying schema..."
psql -h localhost -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -f "$SCRIPT_DIR/../db/schema.sql"

echo "Done. Tables:"
psql -h localhost -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -c "\dt"
