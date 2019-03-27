#!/bin/zsh

# I cycle through all my test images, and make sure they all produce results.
# The grid-finder is relatively anal, so if we found ANYTHING, I assume its good

imagedir=${0:h}/data

failures=()

for image ($imagedir/*)
{
    ./mrgingham --clahe $image > /dev/null || failures+=$image
}

if (($#failures)); then
    # print each failure on a separate line. No idea why this works
    echo "Got $#failures failures:"
    printf " %s\n" "${failures[@]}"
    exit 1
fi

echo "All passed!"
exit 0

   

