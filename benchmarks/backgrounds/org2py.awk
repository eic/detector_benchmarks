BEGIN {
	in_src = 0
	IGNORECASE=1
}

/^#\+begin_src\s+[^\s]*python/ {
	in_src = 1
	match($0, /^ */)
	spaces = RLENGTH
	print "import time; import sys; sys.stderr.write(f\"{time.asctime()}\\n\")"
	next
}

/^#\+end_src/ {
	in_src = 0
	next
}

in_src {
	re = "^ {" spaces "}"
	gsub(re,"")
	print
}
