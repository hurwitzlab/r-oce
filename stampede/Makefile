APP = r-oce-0.0.1
EMAIL = kyclark@email.arizona.edu

clean:
	find . \( -name \*.out -o -name \*.param \) -exec rm {} \;

up: 
	copy_from_manifest.py

app:
	apps-addupdate -F app.json

template:
	jobs-template -A $(APP) > job.json

job:
	jobs-submit -F job.json
