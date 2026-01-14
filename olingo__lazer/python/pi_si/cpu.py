import psutil
import multiprocessing
import time

def cpu_intensive_task():
    """Single-threaded CPU load"""
    while True:
        sum([i**2 for i in range(10000)])

def monitor_single_core(duration=10):
    """Monitor CPU frequency with load on 1 core"""
    print("Starting single-core load test...")
    print(f"Total cores: {psutil.cpu_count()}")
    print(f"Initial frequency: {psutil.cpu_freq().current} MHz\n")
    
    # Start load on ONE core only
    process = multiprocessing.Process(target=cpu_intensive_task)
    process.start()
    
    frequencies = []
    
    try:
        start_time = time.time()
        while time.time() - start_time < duration:
            # Overall frequency
            freq = psutil.cpu_freq()
            
            # Per-core frequency (if available)
            per_core = psutil.cpu_freq(percpu=True)
            
            # CPU usage
            cpu_percent = psutil.cpu_percent(interval=0.1, percpu=True)
            
            elapsed = time.time() - start_time
            print(f"Time: {elapsed:5.1f}s | Overall Freq: {freq.current:7.2f} MHz")
            
            # Show per-core info
            if per_core:
                for i, (core_freq, core_usage) in enumerate(zip(per_core, cpu_percent)):
                    marker = " <-- LOADED" if core_usage > 80 else ""
                    print(f"  Core {i}: {core_freq.current:7.2f} MHz | Usage: {core_usage:5.1f}%{marker}")
            else:
                # If per-core freq not available, just show usage
                for i, core_usage in enumerate(cpu_percent):
                    marker = " <-- LOADED" if core_usage > 80 else ""
                    print(f"  Core {i} Usage: {core_usage:5.1f}%{marker}")
            
            print()
            frequencies.append(freq.current)
            time.sleep(1)
    
    finally:
        print("Stopping load...")
        process.terminate()
        process.join()
        
        # Summary
        if frequencies:
            print(f"\n--- Summary ---")
            print(f"Average frequency: {sum(frequencies)/len(frequencies):.2f} MHz")
            print(f"Peak frequency: {max(frequencies):.2f} MHz")

if __name__ == "__main__":
    monitor_single_core(duration=10)
