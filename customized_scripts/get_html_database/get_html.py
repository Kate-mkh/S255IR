from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from bs4 import BeautifulSoup
import pandas as pd
import time

def get_quantum_numbers(frequency):

    frequency_str = str(frequency)
    
    # setting Selenium and Firefox
    driver = webdriver.Firefox()  # requiered geckodriver
    url = 'https://pml.nist.gov/cgi-bin/micro/table5/start.pl'
    driver.get(url)

    try:
        wait = WebDriverWait(driver, 10)
        input_field = wait.until(EC.presence_of_element_located((By.NAME, 'frequency_area')))
        input_field.clear()
        input_field.send_keys(frequency_str)
        
        # input to 'frequency_range'
        range_field = driver.find_element(By.NAME, 'frequency_range')
        range_field.clear()
        range_field.send_keys("1")
        range_field.send_keys(Keys.RETURN) # Enter
        time.sleep(3) # page loading... 

        # getting HTML 
        html = driver.page_source
        soup = BeautifulSoup(html, 'html.parser')
        table = soup.find('table')
        if not table:
            return "No results"

        rows = table.find_all('tr')
        
        results = []
        for row in rows[1:]:  # skip header
            cols = row.find_all('td')
            if len(cols) > 4:
                results.append(cols[4].text.strip())  # we need line parameters
        
        if len(results) == 1:
            return results[0]
        elif results:
            return " or ".join(results)
        else:
            return "No valid results found"

    finally:
        # Bye browser
        driver.quit()


frequencies = [
    244789.253,
    244832.183,
    243120.317,
    234255.27,
    234758.793,
    246105.97,
    236717.2,
    240982.77,
    240982.77,
    241946.83,
    238156.319,
    216112.58
    ]

# Quantum Numbers for all the freqs
quantum_numbers = [get_quantum_numbers(freq) for freq in frequencies]

df = pd.DataFrame({
    'Frequency (MHz)': frequencies,
    'Quantum Numbers': quantum_numbers
})

df.to_csv('get_html.csv')

print(df['Quantum Numbers'])
