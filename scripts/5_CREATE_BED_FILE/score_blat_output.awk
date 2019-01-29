#!/bin/awk -f
BEGIN {
	OFS="\t"
}
{
	qSize = $13 - $12
	tSize = $17 - $16
	totalSize = qSize > tSize ? qSize : tSize
	score = $1 / totalSize
	if (score > .98) {
		print $14,$16,$17,$10,score
	}
}
