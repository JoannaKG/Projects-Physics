import requests
from bs4 import BeautifulSoup
import re

r = requests.get('http://osn.pwr.edu.pl/group-members.php').content
soup = BeautifulSoup(r, 'lxml')
with open('OSN_nazwiska.csv', 'w') as f:
    for div in soup.find_all("div", {'class': "col-xs-12 col-sm-8 col-md-8"}):
        name = re.sub('Prof. ', '', div.find("h4").text)
        name = re.sub('Dr. ', '', name)
        names = re.findall('(\w+\-*\w+)', name)
        given_name = names[0]
        family_name = names[1]
        mail = div.find("a").text
        f.write(given_name + ', ' + family_name + ', ' + name + ',, ' + mail + ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n')
        print(given_name, family_name, mail)
