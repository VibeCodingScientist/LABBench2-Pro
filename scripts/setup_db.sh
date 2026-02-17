#!/usr/bin/env bash
set -euo pipefail

DB_NAME="${1:-labbench2pro}"
DB_USER="${2:-dev}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Creating database '$DB_NAME' (if not exists)..."
createdb -U "$DB_USER" "$DB_NAME" 2>/dev/null || echo "Database already exists."

echo "Applying schema..."
psql -U "$DB_USER" -d "$DB_NAME" -f "$SCRIPT_DIR/../db/schema.sql"

echo "Done. Tables:"
psql -U "$DB_USER" -d "$DB_NAME" -c "\dt"
