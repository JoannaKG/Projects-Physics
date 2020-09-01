import json


def php_writer(filename, articles):
    """
    Function converts input list of dictionaries with data about articles into .php format to be presented on website

    :param filename: filename of output data
    :param articles: loaded json file
    :return: output file '[filename].php'
    """
    with open(filename, 'w', encoding='UTF-8') as file:
        file.write("<?\necho('\n")
        for article in articles:
            number = '(' + article['issue'] + ')' if 'issue' in article.keys() else ''
            pages = ', ' + article['page'] if 'page' in article.keys() else ''

            file.write('\n<a href="https://doi.org/' + article['DOI'] + '" target="_blank">\n'
                                                                        '<b>' + article['title'] + '</b>\n<br/>')
            file.write(', '.join([author['given'][0] + '. ' + author['family'] for author in article['author']]))
            file.write('<br/>\n<i>' + article['container-title'] + ' ' + article['volume'])
            file.write(number + pages)
            file.write(' (' + article['issued']['date-parts'][0][0] + ')' + '</i>\n')
            file.write('</a>\n\n<br/><br/>\n')
        file.write("\n\n');\n?>")


if __name__ == "__main__":
    with open('2020.json', encoding='UTF-8') as f:
        articles = json.load(f)

    output_file = '2020.php'
    php_writer(output_file, articles)
