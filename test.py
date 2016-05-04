import re
import unicodedata
def slugify(value):
    """
    Normalizes string, converts to lowercase, removes non-alpha characters,
    and converts spaces to hyphens.
    """
    value = unicode(re.sub('[^\w\s-]', '', value).strip().lower())
    value = re.sub('[-\s]+', '-', value)
    return(value)

print(slugify('Phpat.005G030300.1.p pacid=28243865 transcript=Phpat.005G030300.1 locus=Phpat.005G030300 ID=Phpat.005G030300.1.v3.0 annot-version=v3.0'))
