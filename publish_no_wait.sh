#
echo did you refresh the pages by running?
echo 'julia --color=yes make.jl && cp -a build/* site/'
echo 'firefox site/index.html'
echo 'Publishing'
export DOCSDIR=`pwd`
export newdir=/tmp/neppack$$
mkdir $newdir
cd $newdir
git clone -b "gh-pages" git@github.com:nep-pack/NonlinearEigenproblems.jl.git
cd NonlinearEigenproblems.jl/
cp -r $DOCSDIR/site/* .
git add *;  git commit . -m "refresh docs"; git push
cd $DOCSDIR
rm -Rf $newdir
