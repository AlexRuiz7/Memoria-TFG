var express = require('express');
var app = express();

//MESSAGES ENDPOINT
app.post('/deviceMessage', function(req, res) {

    // Extract and separate data and cryptographic signature
    var crypto = require('crypto');
    var verify = crypto.createVerify('RSA-SHA256');
    var d = JSON.parse(JSON.stringify(req.body));
    var sign = new Buffer(d.sign, 'base64');
    delete d.sign
    verify.write(d.data);
    verify.end();
    var message = JSON.parse(d.data)

    //Searching for the device that is trying to save data
    db.get("SELECT * FROM devices WHERE serialno = ?", message.serialno, function(err, row) {

        //If device exists and is registered
        if (row && row.public_key) {

            //Verify the cryptographic signature
            if (verify.verify(row.public_key, sign)) {
                //Insert the data in the database
                db.run("-SQL DATA INSERTION-", [data], function(err) {
                    console.log("HTTP Measure saved as " + this.lastID)
                    res.json({status: "OK"})
                });
            //Not valid cryptographic signature
            } else {
                res.json({error: "Signature Error"})
            }
        //Nonexistent or not registered device
        } else {
            res.json({error: "Non valid serialno"})
        }
    })
});

//REGISTRATION ENDPOINT
app.post('/deviceRegister', function(req, res) {
    var data = req.body

    //Search the device that is trying to save data in the system
    db.get("SELECT * FROM devices WHERE serialno = ?", data.serialno, function(err, row) {

        //If device exists
        if (row) {
            //If device is not registered
            if (row.public_key == null) {
                db.run("-SQL PUBLIC_KEY INSERTION-",data.public_key, function(err) {
                    res.json({status: "OK"})
                })
            //Device already registered
            } else {
                res.json({error: "Registered device"})
            }
        //Nonexistent or not registered device
        } else {
            res.json({error: "Non valid serial number"})
        }
    })
});

//Start the webserver in a certain TCP port
app.listen(WEBPORT, HOST);
