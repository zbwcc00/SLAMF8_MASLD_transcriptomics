#!/usr/bin/env python
import re
import requests


def main():
    for accession in ["GSE135251", "GSE126848", "GSE162694", "GSE130970", "GSE193066", "GSE167523"]:
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/suppl/"
        response = requests.get(url, timeout=35)
        names = re.findall(r'href="([^"]+)"', response.text)
        print(accession, response.status_code, names[-15:])


if __name__ == "__main__":
    main()
