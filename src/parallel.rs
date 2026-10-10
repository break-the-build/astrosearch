use crate::Candidate;
use std::{
    ops::Range,
    sync::atomic::{AtomicUsize, Ordering},
    thread,
};

const CHUNK_SIZE: usize = 16;

/// Independent trials retain their serial ordering before scientific ranking.
/// Dynamic chunks balance slow/fast cores and the varying cost of trial periods.
/// Each worker owns and reuses its state; light-curve arrays remain shared.
pub(crate) fn collect_trials<State: Send>(
    trials: usize,
    requested_threads: usize,
    initialize: impl Fn() -> State + Sync,
    evaluate: impl Fn(&mut State, usize) -> Vec<Candidate> + Sync,
) -> Vec<Candidate> {
    let available = if requested_threads == 0 {
        thread::available_parallelism().map_or(1, usize::from)
    } else {
        requested_threads
    };
    let workers = available.min(trials.div_ceil(CHUNK_SIZE)).max(1);
    if workers == 1 {
        let mut state = initialize();
        return (0..trials)
            .flat_map(|index| evaluate(&mut state, index))
            .collect();
    }
    let next = AtomicUsize::new(0);
    let mut indexed = thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            let (initialize, evaluate, next) = (&initialize, &evaluate, &next);
            handles.push(scope.spawn(move || {
                let mut state = initialize();
                let mut output = Vec::new();
                loop {
                    let start = next.fetch_add(CHUNK_SIZE, Ordering::Relaxed);
                    if start >= trials {
                        break;
                    }
                    for index in start..(start + CHUNK_SIZE).min(trials) {
                        for (ordinal, candidate) in
                            evaluate(&mut state, index).into_iter().enumerate()
                        {
                            output.push((index, ordinal, candidate));
                        }
                    }
                }
                output
            }));
        }
        let mut output = Vec::new();
        for handle in handles {
            output.extend(handle.join().expect("search worker panicked"));
        }
        output
    });
    indexed.sort_unstable_by_key(|(index, ordinal, _)| (*index, *ordinal));
    indexed
        .into_iter()
        .map(|(_, _, candidate)| candidate)
        .collect()
}

/// Process contiguous trial blocks so algorithms with a recurrence across
/// adjacent frequencies can reuse phase state while retaining parallel workers.
pub(crate) fn collect_trial_blocks<State: Send>(
    trials: usize,
    block_size: usize,
    requested_threads: usize,
    initialize: impl Fn() -> State + Sync,
    evaluate: impl Fn(&mut State, Range<usize>) -> Vec<(usize, Candidate)> + Sync,
) -> Vec<Candidate> {
    let block_size = block_size.max(1);
    let blocks = trials.div_ceil(block_size);
    let available = if requested_threads == 0 {
        thread::available_parallelism().map_or(1, usize::from)
    } else {
        requested_threads
    };
    let workers = available.min(blocks).max(1);
    let mut indexed = if workers == 1 {
        let mut state = initialize();
        (0..blocks)
            .flat_map(|block| {
                let start = block * block_size;
                evaluate(&mut state, start..(start + block_size).min(trials))
            })
            .collect()
    } else {
        let next = AtomicUsize::new(0);
        thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for _ in 0..workers {
                let (initialize, evaluate, next) = (&initialize, &evaluate, &next);
                handles.push(scope.spawn(move || {
                    let mut state = initialize();
                    let mut output = Vec::new();
                    loop {
                        let block = next.fetch_add(1, Ordering::Relaxed);
                        if block >= blocks {
                            break;
                        }
                        let start = block * block_size;
                        output.extend(evaluate(
                            &mut state,
                            start..(start + block_size).min(trials),
                        ));
                    }
                    output
                }));
            }
            let mut output = Vec::new();
            for handle in handles {
                output.extend(handle.join().expect("search worker panicked"));
            }
            output
        })
    };
    indexed.sort_unstable_by_key(|(index, _)| *index);
    indexed
        .into_iter()
        .map(|(_, candidate)| candidate)
        .collect()
}
