#!/bin/sh

# Drives a series of simulation experiments in the big_experiments subdirectory
# exploring different settings for the number of tandem readas to generate
# for a range of numbers of input reads

python explore_variations.py \
    --name train_series \
    --start-from beginning == \
    --trials=10 --predict-for-training --keep-intermediates \
    --assess-accuracy --skip-rewrite --sim-unp-min=1 \
    --sim-conc-min=1 --sim-disc-min=1 --sim-bad-end-min=1 \
    == ts_01m_15s  ts_01m_30s  ts_01m_45s \
       ts_01m_1l   ts_01m_3l   ts_01m_5l \
       ts_01m_050c ts_01m_100c \
       ts_05m_15s  ts_05m_30s  ts_05m_45s \
       ts_05m_1l   ts_05m_3l   ts_05m_5l \
       ts_05m_050c ts_05m_100c \
       ts_10m_15s  ts_10m_30s  ts_10m_45s \
       ts_10m_1l   ts_10m_3l   ts_10m_5l \
       ts_10m_050c ts_10m_100c \
       ts_50m_15s  ts_50m_30s  ts_50m_45s \
       ts_50m_1l   ts_50m_3l   ts_50m_5l \
       ts_50m_050c ts_50m_100c \
    == --sim-function sqrt --sim-factor 15 -- --sim-function sqrt --sim-factor 30 -- --sim-function sqrt --sim-factor 45 \
    -- --sim-function linear --sim-factor 0.01 -- --sim-function linear --sim-factor 0.03 -- --sim-function linear --sim-factor 0.05 \
    -- --sim-function const --sim-factor 50000 -- --sim-function const --sim-factor 100000 \
    -- --sim-function sqrt --sim-factor 15 -- --sim-function sqrt --sim-factor 30 -- --sim-function sqrt --sim-factor 45 \
    -- --sim-function linear --sim-factor 0.01 -- --sim-function linear --sim-factor 0.03 -- --sim-function linear --sim-factor 0.05 \
    -- --sim-function const --sim-factor 50000 -- --sim-function const --sim-factor 100000 \
    -- --sim-function sqrt --sim-factor 15 -- --sim-function sqrt --sim-factor 30 -- --sim-function sqrt --sim-factor 45 \
    -- --sim-function linear --sim-factor 0.01 -- --sim-function linear --sim-factor 0.03 -- --sim-function linear --sim-factor 0.05 \
    -- --sim-function const --sim-factor 50000 -- --sim-function const --sim-factor 100000 \
    -- --sim-function sqrt --sim-factor 15 -- --sim-function sqrt --sim-factor 30 -- --sim-function sqrt --sim-factor 45 -- --sim-function linear --sim-factor 0.01 -- --sim-function linear --sim-factor 0.03 -- --sim-function linear --sim-factor 0.05 -- --sim-function const --sim-factor 50000 -- --sim-function const --sim-factor 100000 \
    == -u 1000000 -- -u 1000000 -- -u 1000000 \
    -- -u 1000000 -- -u 1000000 -- -u 1000000 \
    -- -u 1000000 -- -u 1000000 \
    -- -u 5000000 -- -u 5000000 -- -u 5000000 \
    -- -u 5000000 -- -u 5000000 -- -u 5000000 \
    -- -u 5000000 -- -u 5000000 \
    -- -u 10000000 -- -u 10000000 -- -u 10000000 \
    -- -u 10000000 -- -u 10000000 -- -u 10000000 \
    -- -u 10000000 -- -u 10000000 \
    -- -u 50000000 -- -u 50000000 -- -u 50000000 \
    -- -u 50000000 -- -u 50000000 -- -u 50000000 \
    -- -u 50000000 -- -u 50000000 \
    == big_experiments/r0_bt2s_mason_ill_100_50M.out big_experiments/r12_bt2s100_mason_ill_100_50M.out
