from pattern.web import *

g = Google()
for result in g.search('lucy'):
	print result.url
	print plaintext(result.text)