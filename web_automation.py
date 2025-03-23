# -*- coding: utf-8 -*-
import os
import time
import random
import glob
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

def setup_driver(download_dir):
    chrome_options = webdriver.ChromeOptions()
    prefs = {
        "download.default_directory": download_dir,
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing.enabled": True
    }
    chrome_options.add_experimental_option("prefs", prefs)
    chrome_options.add_argument("--disable-blink-features=AutomationControlled")
    return webdriver.Chrome(
        service=Service(ChromeDriverManager().install()),
        options=chrome_options
    )

def process_sequence(driver, download_dir, sequence, seq_id, dna_type):
    try:
        driver.get("https://mscskeylab.hus.vnu.edu.vn/biotools/dnagenerator.php")
        
        # Fill form
        WebDriverWait(driver, 15).until(
            EC.presence_of_element_located((By.NAME, "name"))
        ).send_keys(seq_id[:4])

        Select(driver.find_element(By.NAME, "dnatype")).select_by_visible_text(
            "Right Handed A-DNA (Arnott)" if dna_type == "A" 
            else "Right Handed B-DNA (Arnott)"
        )

        seq_field = driver.find_element(By.NAME, "DNAsequence")
        seq_field.clear()
        seq_field.send_keys(sequence.strip())

        # Submit form
        driver.find_element(By.XPATH, "//input[@type='submit']").click()

        # Handle download
        WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.XPATH, "//a[contains(text(),'Download')]"))
        ).click()

        # Wait for download completion
        max_wait = 30
        start_time = time.time()
        while (time.time() - start_time) < max_wait:
            files = glob.glob(os.path.join(download_dir, "*.pdb"))
            if files and os.path.getsize(max(files, key=os.path.getctime)) > 0:
                break
            time.sleep(1)

        # Rename downloaded file
        latest_file = max(files, key=os.path.getctime)
        new_name = os.path.join(download_dir, f"{seq_id}_{dna_type}.pdb")
        os.rename(latest_file, new_name)

        return True

    except Exception as e:
        print(f"Error processing {seq_id}: {str(e)}")
        driver.save_screenshot(f"error_{seq_id}_{dna_type}.png")
        return False
    finally:
        # Reset to initial page
        driver.get("https://mscskeylab.hus.vnu.edu.vn/biotools/dnagenerator.php")
        time.sleep(random.uniform(1, 2))

def dna_sequence_processor(csv_path, a_dna_dir, b_dna_dir):
    df = pd.read_csv(csv_path)
    
    # Create separate drivers for each DNA type
    a_driver = setup_driver(a_dna_dir)
    b_driver = setup_driver(b_dna_dir)

    try:
        for index, row in df.iterrows():
            seq = row['Sequence']
            seq_id = f"seq{index:04d}"

            # Process A-DNA
            process_sequence(a_driver, a_dna_dir, seq, seq_id, "A")
            
            # Process B-DNA
            process_sequence(b_driver, b_dna_dir, seq, seq_id, "B")
            
            # Random delay between sequences
            time.sleep(random.uniform(2, 4))

    finally:
        a_driver.quit()
        b_driver.quit()

if __name__ == "__main__":
    # Example usage
    dna_sequence_processor(
        "random_sequence.csv",
        "A_DNA_Structures",
        "B_DNA_Structures"
    )