import re


def bibtex_reader(data, parameters):
    """
    Function parse input data in .bibtex format into dictionary

    :param data: input data in .bibtex format
    :param parameters: dictionary for which new key : value pairs will be added based on parsed data;
                       keys are equal to .bibtex keys and values contain parsed values for the keys from data
    :return:
    is_new_article  - boolean; True if parsed keyword is '@article' (means that we start parsing data about next article)
    """
    is_new_article = False
    if re.findall('(@article)', data):
        article = re.findall('@\w+{(\w+)', data)[0]
        parameters['article'] = article
        is_new_article = True
    elif re.findall('(\t)', data):
        key = re.findall('\t(\w+)', data)[0]
        if key == 'title':
            title = re.findall('\ttitle = (.*),', data)[0]
            title = re.sub('{|}|\\\\textbf|\\\\textit|\\\\', '', title)
            parameters['title'] = title
        elif key == 'author':
            author = re.findall('\tauthor = (.*),', data)[0]
            authors = []
            author = re.findall('\{*(\w*-*\w*?), (\w+)', author)
            for surname, name in author:
                authors.append(name[0] + '. ' + surname)
            author = ', '.join(authors)
            parameters['author'] = author
        elif key == 'month':
            pass
        else:
            value = re.findall('\t' + key + ' = {(.*)},*', data)[0]
            value = re.sub('--', '-', value)
            parameters[key] = value
    else:
        pass
    return is_new_article


def php_writer(filename, articles):
    """
    Function converts input list of dictionaries with data about articles into .php format to be presented on website

    :param filename: filename of output data
    :param articles: list of dictionaries containing keys equal to .bibtex keys and values containing parsed values
    :return: output file '[filename].php'
    """
    with open(filename, 'w', encoding='UTF-8') as file:
        file.write("<?\necho('\n")
        for article in articles:
            number = '(' + article['number'] + ')' if 'number' in article.keys() else ''
            pages = ', ' + article['pages'] if 'pages' in article.keys() else ''

            file.write('\n<a href="https://doi.org/' + article['doi'] + '" target="_blank">\n'
                                                                        '<b>' + article['title'] + '</b>\n<br/>' +
                       article['author'] + '<br/>')
            file.write('\n<i>' + article['journal'] + ' ' + article['volume'])
            file.write(number + pages)
            file.write(' (' + article['year'] + ')' + '</i>\n')
            file.write('</a>\n\n<br/><br/>\n')
        file.write("\n\n');\n?>")


if __name__ == "__main__":
    input_file = '2020.bib'
    articles = []
    parameters = {}
    with open(input_file, encoding='UTF-8') as file:
        data = file.readlines()
        for line in data:
            if bibtex_reader(line, parameters):
                parameters = {}
                articles.append(parameters)
    output_file = '2020.php'
    php_writer(output_file, articles)
