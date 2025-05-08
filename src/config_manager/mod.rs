//
// config_manager/main.rs: Handle data files and defult paths
//
// Copyright (c) 2025 Pablo Álvarez Martín <pablo.alvmar12@gmail.com>
// Copyright (c) 2025 Gonzalo J. Carracedo <BatchDrake@gmail.com>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//

use crate::data_file_manager::*;
use serde_yaml::Value; // this is for the YAML nodes as a generic type
use std::collections::HashMap;
use std::fs; // for create_dir (equivalent of mkdir)

use nix::errno::Errno; 
use std::ffi::CString; 
use nix::libc::*;
use std::path::Path;
use errno::errno;
use once_cell::sync::Lazy; 
use std::sync::Mutex; // for mutable instance

use std::fs::File;
use std::io::Read;
use std::io::Write;

// use lazy_static::lazy_static;

//use serde::de::DeserializeOwned;
use serde::Serialize;

// lazy_static! {
//     static ref CONFIG_MANAGER_DIRECTORY: String = String::from("/ZarzaETC/config");
// }

#[derive(Default, Debug, Clone, PartialEq)]
pub struct Config {
    pub object_name: String,
    pub yaml_config: Value,
}

impl Config {

    //pub fn config(&mut self, name: String) {
    pub fn new(name: String) -> Self {

        let mut config = Config {
            object_name: String::new(),
            yaml_config: Value::Null,
        };

        if name.is_empty() {
            eprintln!("Empty config names are not allowed.");
        }
        config.object_name = name;

        let _p= &config.yaml_config["a"]; // ???

        return config;
    }

    pub fn serialize() -> bool { // ???????
        return true;
    }

    pub fn deserialize() -> bool { // ???????
        return true;
    }

    pub fn load(&mut self) -> bool {
        let manager = ConfigManager::instance().lock().unwrap();
        let path = manager.get_config_file_path(self.object_name.clone(), false);

        if Path::new(&path).exists() { // Check if file exists
            let mut file = File::open(path).expect("File cannot be opened."); // Check if file can be read (and open it)
            let mut yaml_content = String::new();
            file.read_to_string(&mut yaml_content).expect("Failed to read file");
            self.yaml_config = serde_yaml::from_str(&yaml_content).expect("Failed to parse YAML");
        } else if !Path::new(&path).exists() {
            eprintln!("{} does not exist.", path);
        };
        
        return Config::deserialize();
    }

    pub fn save(&self) -> bool {
        let manager = ConfigManager::instance().lock().unwrap();
        let path = manager.get_config_file_path(self.object_name.clone(), false);
    
        let mut file = File::create(&path).expect("Failed to create file.");
        let yaml_string = serde_yaml::to_string(&self.yaml_config).expect("Failed to write string");
        file.write_all(yaml_string.as_bytes()).expect("Fail to write string into .yaml");

        return Config::serialize();
    }

    pub fn has_key(&self, name: String) -> bool {
        return self.yaml_config.get(name).is_some();
    }

    pub fn yaml_node(&self) -> Value {
        return self.yaml_config.clone();
    }

    pub fn deserialize_yaml_node(&mut self, node: &Value) -> bool {
        self.yaml_config = node.clone();
        return Config::deserialize();
    }

    pub fn serialize_field<T>(field: &T) -> Value // substitution of "operator[]"
    where 
        T: Serialize 
    {
        serde_yaml::to_value(field).expect("Failed to serialize field")
    }

    pub fn deserialize_field<T>(&mut self, dest: &mut T, name: &str) -> bool 
    where
        T: serde::de::DeserializeOwned, // This adds the required bound
    {

        let value = match self.yaml_config.get(&Value::String(name.to_string())) {
            Some(v) => v,
            None => return false, // Okay: use default (like C++ logic)
        };
        // deserialize field into the target type:
        let val = serde_yaml::from_value(value.clone()).expect("Failed to deserialize key."); // Check error handling with expects, find alternatives or is ok for the program to panic?
        *dest = val;
        return true;

    }


}

pub struct ConfigManager {
    config_dir: String,
    config_list: Vec<Config>,
    config_cache: HashMap<String, Config>,
    can_save_config: bool,
}

impl ConfigManager {

    //pub fn config_manager(&mut self) {
    pub fn new() -> Self {

        let mut config_manager = ConfigManager {
            config_dir: String::new(),
            config_list: Vec::new(),
            config_cache: HashMap::new(),
            can_save_config: false, //??
        };

        let mut sbuf: stat = unsafe { std::mem::zeroed() }; 

        //config_manager.config_dir = DataFileManager::instance().suggest(CONFIG_MANAGER_DIRECTORY.clone());
        config_manager.config_dir = DataFileManager::instance().suggest("config".to_string());

        if !config_manager.config_dir.is_empty() {
            let c_config_dir = CString::new(config_manager.config_dir.clone()).expect("CString::new failed");
            if unsafe { stat(c_config_dir.as_ptr(), &mut sbuf as *mut stat) } == -1 { // check "metadata", might be safer than stat...
                if Errno::last() == Errno::ENOENT {
                    let path = Path::new(&config_manager.config_dir);
                    match fs::create_dir(path) {
                        Ok(_) => {
                            config_manager.can_save_config = true;
                        }
                        Err(e) => {
                            eprintln!("Failed to create config directory: {}", e);
                        }  
                    }                    
                } else {
                    let err = errno();
                    eprintln!("directory `{}`inaccessible: {}", config_manager.config_dir, err.to_string());
                }
            } else if sbuf.st_mode & S_IFDIR == 0 {
                let err = errno();
                eprintln!("directory `{}`inaccessible: {}", config_manager.config_dir, err.to_string());
            // I think the "catch" part of the c++ code is not necessary because of rust error handling (?)
            } else {
                config_manager.can_save_config = true;
            }
        }
        if config_manager.can_save_config == false {
            eprintln!("warning: no writable config directory available, configurations cannot be saved!"); 
        }

        return config_manager;
    }

    pub fn save_all(&self) -> bool {

        let mut ok: bool = true;

        if !self.can_save_config {
            return false
        }

        for config in &self.config_list {
            ok = config.save() && ok;
        }

        return ok;
    }

    pub fn get_config_file_path(&self, name: String, write: bool) -> String {
        if write == true {
            return self.config_dir.clone() + "/" + &name + ".yaml";
        } else {
            //return DataFileManager::instance().resolve(CONFIG_MANAGER_DIRECTORY.clone() + "/" + &name + ".yaml");
            return DataFileManager::instance().resolve("config".to_owned() + "/" + &name + ".yaml");
        }
    }

    pub fn instance() -> &'static Mutex<ConfigManager> { // Necessary to use Mutex to make the instance mutable
        static INSTANCE: Lazy<Mutex<ConfigManager>> = Lazy::new(|| Mutex::new(ConfigManager::new()));
    
        // Return a reference to the Mutex, which can be locked for mutable access
        &INSTANCE
    }

    // I dont use generic types because inheritance doesnt work the same as in c++, 
    // instead i will implement these functions in Detector, Simulation, ... as traits of the 
    // objects derived from Config and ConfigManager that are in detector, simulation, etc. 

    pub fn get_config(&mut self, name: String) -> &Config { 

        if !self.config_cache.contains_key(&name) {
            let mut new_config = Config::new(name.clone());
            new_config.load();
            self.config_list.push(new_config.clone());
            self.config_cache.insert(name.clone(), new_config);
        }
        // no else condition necessary, implies that Config was found, and is returned below.

        return self.config_cache.get(&name).unwrap();
    }

    pub fn get(name: String) -> Config {
        let mut manager = ConfigManager::instance().lock().unwrap();

        let cfg = manager.get_config(name);

        if cfg.object_name.is_empty() {
            eprintln!("Failed to retrieve config.");
        }

        return cfg.clone();
    }
} 



