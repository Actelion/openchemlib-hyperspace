package com.idorsia.research.chem.hyperspace.util;

import java.util.List;
import java.util.concurrent.*;
import java.util.function.Consumer;

public class Parallelizer {

    public static <T> void computeParallelBlocking(Consumer<T> f, List<T> objects, int threads) throws InterruptedException {
        ExecutorService executor = Executors.newFixedThreadPool(threads);

        for (T obj : objects) {
            executor.submit(() -> f.accept(obj));
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

    static class Stats {
        private int completedTasks;
        private long startTime;
        private long lastOutputTime;

        public Stats(int completedTasks, long startTime, long lastOutputTime) {
            this.completedTasks = completedTasks;
            this.startTime = startTime;
            this.lastOutputTime = lastOutputTime;
        }

        public synchronized int getCompletedTasks() {
            return completedTasks;
        }

        public synchronized void incrementCompletedTasks() {
            completedTasks++;
        }

        public synchronized long getStartTime() {
            return startTime;
        }

        public synchronized long getLastOutputTime() {
            return lastOutputTime;
        }

        public synchronized void setLastOutputTime(long lastOutputTime) {
            this.lastOutputTime = lastOutputTime;
        }
    }

    public static void computeParallelBlocking(List<Runnable> tasks, int threads, int outputIntervalSeconds) throws InterruptedException {
        ExecutorService executor = Executors.newFixedThreadPool(threads);
        CompletionService<Void> completionService = new ExecutorCompletionService<>(executor);

        int totalTasks = tasks.size();
        Stats stats = new Stats(0, System.currentTimeMillis(), System.currentTimeMillis());

        for (Runnable task : tasks) {
            completionService.submit(task, null);
        }

        executor.shutdown();

        ScheduledExecutorService outputScheduler = Executors.newScheduledThreadPool(1);
        outputScheduler.scheduleAtFixedRate(() -> {
            long currentTime = System.currentTimeMillis();
            long elapsedTime = currentTime - stats.getStartTime();
            double currentSpeed = (double) stats.getCompletedTasks() / elapsedTime;
            long remainingTasks = totalTasks - stats.getCompletedTasks();
            double expectedRemainingTime = remainingTasks / currentSpeed;

            System.out.printf("Completed tasks: %d/%d | Current speed: %.2f tasks/ms | Expected remaining time: %.2f ms\n",
                    stats.getCompletedTasks(), totalTasks, currentSpeed, expectedRemainingTime);

            stats.setLastOutputTime(currentTime);
        }, outputIntervalSeconds * 1000, outputIntervalSeconds * 1000, TimeUnit.MILLISECONDS);

        while (stats.getCompletedTasks() < totalTasks) {
            Future<Void> completedTask = completionService.take();
            stats.incrementCompletedTasks();
        }

        outputScheduler.shutdown();
        executor.shutdown();

        // Print final statistics after all tasks are completed
        long currentTime = System.currentTimeMillis();
        long elapsedTime = currentTime - stats.getStartTime();
        double currentSpeed = (double) stats.getCompletedTasks() / elapsedTime;
        long remainingTasks = totalTasks - stats.getCompletedTasks();
        double expectedRemainingTime = remainingTasks / currentSpeed;

        System.out.printf("Completed tasks: %d/%d | Current speed: %.2f tasks/ms | Expected remaining time: %.2f ms\n",
                stats.getCompletedTasks(), totalTasks, currentSpeed, expectedRemainingTime);
    }

}
