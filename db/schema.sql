CREATE TABLE IF NOT EXISTS models (
    id          SERIAL PRIMARY KEY,
    provider    TEXT NOT NULL,
    model_name  TEXT NOT NULL,
    tools       BOOLEAN DEFAULT FALSE,
    UNIQUE(provider, model_name, tools)
);

CREATE TABLE IF NOT EXISTS tasks (
    id          TEXT PRIMARY KEY,
    category    TEXT NOT NULL,
    variant     TEXT,
    source      TEXT DEFAULT 'labbench2',
    meta        JSONB
);

CREATE TABLE IF NOT EXISTS eval_runs (
    id          SERIAL PRIMARY KEY,
    model_id    INTEGER REFERENCES models(id),
    task_id     TEXT REFERENCES tasks(id),
    response    TEXT,
    correct     BOOLEAN,
    score       FLOAT,
    grader      TEXT,
    tokens_in   INTEGER,
    tokens_out  INTEGER,
    cost_usd    FLOAT,
    latency_ms  INTEGER,
    run_at      TIMESTAMPTZ DEFAULT NOW(),
    meta        JSONB
);

CREATE TABLE IF NOT EXISTS chain_runs (
    id          SERIAL PRIMARY KEY,
    chain_id    TEXT NOT NULL,
    model_id    INTEGER REFERENCES models(id),
    step_num    INTEGER NOT NULL,
    task_id     TEXT REFERENCES tasks(id),
    input_from  TEXT,
    response    TEXT,
    correct     BOOLEAN,
    run_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE IF NOT EXISTS judge_audits (
    id           SERIAL PRIMARY KEY,
    eval_run_id  INTEGER REFERENCES eval_runs(id),
    judge_model  TEXT NOT NULL,
    judge_score  BOOLEAN,
    order_variant TEXT,
    length_variant TEXT,
    explanation  TEXT,
    run_at       TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_eval_model_task ON eval_runs(model_id, task_id);
CREATE INDEX idx_eval_task_correct ON eval_runs(task_id, correct);
CREATE INDEX idx_chain_model ON chain_runs(chain_id, model_id);
CREATE INDEX idx_tasks_category ON tasks(category);
