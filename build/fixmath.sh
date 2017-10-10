mv build/index.md /tmp/index.md
cat /tmp/index.md |  perl -pe 's/([^\$])\$([^\$]+)\$/$1\\\\($2\\\\)/g'  > build/index.md

