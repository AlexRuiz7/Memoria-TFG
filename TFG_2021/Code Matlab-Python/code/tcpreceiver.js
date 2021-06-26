var net = require('net');
var path = require("path");
var dateFormat = require('dateformat');
var pad = require("left-pad")

//Parsing raw data message
function parseMessage(messageData) {
    var data = {}

    var modelCAP = messageData.toString("utf-8", 0, 6)
    var type = messageData.toString("utf-8", 15, 16)

    //Check if CAP21 string is in header and if X is the type
    if (modelCAP == " CAP21" && type == "X") {

        //Extract the serial number
        data.serialno = messageData.toString("utf-8", 7, 15)

        //Extract the message timestamp
        data.date =
        pad(parseInt(messageData.toString("hex",18,19),16),2) + "-" +
        pad(parseInt(messageData.toString("hex",17,18),16),2) + "-" +
        pad(parseInt(messageData.toString("hex",16,17),16),2) + " " +
        pad(parseInt(messageData.toString("hex",20,21),16),2) + ":" +
        pad(parseInt(messageData.toString("hex",19,20),16),2) + ":00"

        //Extract acoustic noise info
        data.L90 = parseInt(messageData.toString("hex",23,24),16)
        data.Leq = parseInt(messageData.toString("hex",24,25),16)
        data.L10 = parseInt(messageData.toString("hex",25,26),16)
        data.Lmax = parseInt(messageData.toString("hex",26,27),16)
    }
    return data
}

//Listening function with callback and config
exports.listen = function(cb, config) {
    var HOST = config.host || '0.0.0.0';
    var PORT = config.port || 12000;
    net.createServer(function(sock) {

        //Device sending data
        sock.on('data', function(messageData) {

            //Save connection info
            var info = {}
            info.address = this.remoteAddress;
            info.port = this.remotePort;
            info.timestamp = dateFormat(new Date(), "yyyy-mm-dd HH:MM:ss")

            //Call the callback with connection info and parsed data
            cb(mergeObjects(info, parseMessage(messageData)))
        });
    }).listen(PORT, HOST);
}
