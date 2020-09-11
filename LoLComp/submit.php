<?php
     
$json = $_POST['jsonData'];
//  turn stringified json object into a PHP object and handle data for fwrite text formatting if necessary
$jsonObj = json_decode($json);

if ($json != NULL && is_string($json)) {
    print("json valid writing to file");
    $file = fopen("../../results/LoLPilotResults.json", "a+");                                                                                           
    fwrite($file, $json);
    fwrite($file, ",");
    fclose($file);

    //  alternative using file_put_contents
    //file_put_contents("results/new_map_data.json", $json, FILE_APPEND | LOCK_EX);
}
else {
    // user has posted invalid JSON, handle the error
    print("json variable is null!");
    $file = fopen("../../results/LoLfailed.txt", "a+");
    foreach ($_POST as $key => $value) {
        $blah = $key. " ". $value;
        fwrite($file, $blah);
    }

    fclose($file);
}
?>
